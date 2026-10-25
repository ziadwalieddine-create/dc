"""Repository policy check: nothing that must stay out of git is tracked, and the visibility is stated.

Checks the files tracked in the current checkout:
  1. size            tracked files above --max-mb (large datasets belong in data/external/, which is git-ignored, or on a data branch)
  2. personal data   raw-sequencing and genotype file types (BAM/CRAM/FASTQ/VCF/BCF, consumer-genotype exports)
  3. credentials     private keys and common token formats (matches are masked in the output)
  4. visibility      repository visibility read live with `gh api` (informational; --require-private makes a public or unknown repository a failure)

Justified exceptions go in inventory/policy_allowlist.json as {"path": "...", "reason": "..."}; each is printed.
Exit status is non-zero if any check fails.

  python3 inventory/check_repo_policy.py [--max-mb 25] [--require-private]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
PERSONAL = re.compile(r'(\.(bam|cram|sam|bai|crai|fastq|fq|bcf|vcf|gvcf)(\.gz)?$)|(\.g\.vcf(\.gz)?$)|(^|/)genome_[^/]*_full_[^/]*\.txt$|(^|/)ancestrydna[^/]*\.txt$|(^|/)23andme[^/]*', re.I)
SECRETS = [
    ('private key block', re.compile(r'-----BEGIN [A-Z ]*PRIVATE KEY-----')),
    ('GitHub token', re.compile(r'\b(ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{30,}\b|\bgithub_pat_[A-Za-z0-9_]{30,}\b')),
    ('AWS access key id', re.compile(r'\bAKIA[0-9A-Z]{16}\b')),
    ('Slack token', re.compile(r'\bxox[baprs]-[A-Za-z0-9-]{10,}\b')),
    ('bearer/JWT token', re.compile(r'\beyJ[A-Za-z0-9_-]{15,}\.[A-Za-z0-9_-]{15,}\.[A-Za-z0-9_-]{10,}\b')),
    ('Google API key', re.compile(r'\bAIza[0-9A-Za-z_-]{35}\b')),
    ('credential assignment', re.compile(r'(?i)\b(password|passwd|secret|api[_-]?key|access[_-]?token)\b\s*[:=]\s*["\'][^"\'\s]{12,}["\']')),
]


def tracked():
    out = subprocess.run(['git', '-C', ROOT, 'ls-files', '-z'], capture_output=True, check=True).stdout.decode()
    return [p for p in out.split('\0') if p]


def visibility():
    if not shutil.which('gh'):
        return None, 'gh CLI not available'
    url = subprocess.run(['git', '-C', ROOT, 'remote', 'get-url', 'origin'], capture_output=True, text=True).stdout.strip()
    m = re.search(r'github\.com[:/]([^/]+)/([^/]+?)(?:\.git)?$', url)
    if not m:
        return None, 'origin is not a github.com remote'
    r = subprocess.run(['gh', 'api', f'repos/{m.group(1)}/{m.group(2)}', '--jq', '.visibility'], capture_output=True, text=True, timeout=60)
    if r.returncode != 0:
        return None, f'gh api failed: {r.stderr.strip()[:120]}'
    return r.stdout.strip(), f'{m.group(1)}/{m.group(2)}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-mb', type=float, default=25)
    ap.add_argument('--require-private', action='store_true')
    a = ap.parse_args()
    ap_path = os.path.join(HERE, 'policy_allowlist.json')
    allow = {x['path']: x['reason'] for x in json.load(open(ap_path))} if os.path.exists(ap_path) else {}
    files = tracked()
    problems = []
    for p in files:
        fp = os.path.join(ROOT, p)
        if not os.path.isfile(fp):
            continue
        size = os.path.getsize(fp)
        if p in allow:
            print(f'allowlisted: {p} ({allow[p]})')
            continue
        if size > a.max_mb * 1e6:
            problems.append(f'size: {p} is {size / 1e6:.1f} MB (limit {a.max_mb} MB)')
        if PERSONAL.search(p):
            problems.append(f'personal-genome file type: {p}')
        if size <= 2_000_000:
            try:
                text = open(fp, encoding='utf-8', errors='strict').read()
            except (UnicodeDecodeError, OSError):
                continue
            for name, rx in SECRETS:
                for m in rx.finditer(text):
                    line = text.count('\n', 0, m.start()) + 1
                    shown = m.group(0)[:4] + '***'
                    problems.append(f'credential pattern ({name}) in {p}:{line} -> {shown}')
    vis, detail = visibility()
    print(f'{len(files)} tracked files checked')
    print(f'visibility: {vis or "unknown"} ({detail})')
    if a.require_private and vis != 'private':
        problems.append(f'repository visibility is {vis or "unknown"}, but --require-private was given')
    if vis == 'public':
        print('NOTE: the repository is public. Everything tracked is world-readable, including history; keep restricted or personal material out.')
    for p in problems:
        print('PROBLEM', p)
    print('policy check:', 'FAIL' if problems else 'PASS')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
