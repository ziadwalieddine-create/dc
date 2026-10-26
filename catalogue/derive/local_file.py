"""Store a local file (for example a licence read from a sanctioned read-only clone) as evidence, with its provenance.

  python3 catalogue/derive/local_file.py SOURCE NAME PATH --url URL [--note TEXT] [--lines A:B]

Writes catalogue/evidence/SOURCE/NAME.json with the text (all of it, or the 1-based inclusive line range A:B), the SHA-256 of the whole file, and, if the
file lies in a git work tree, the repository's origin URL and commit. A source record then cites it as "SOURCE/NAME". Use it only where the file cannot be
fetched over HTTP by catalogue/probe.py (for example repositories that are available to the session as read-only clones).
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def git(path, *args):
    try:
        return subprocess.run(['git', '-C', os.path.dirname(os.path.abspath(path)), *args], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source')
    ap.add_argument('name')
    ap.add_argument('path')
    ap.add_argument('--url', required=True)
    ap.add_argument('--note', default='')
    ap.add_argument('--lines')
    a = ap.parse_args()
    raw = open(a.path, 'rb').read()
    lines = raw.decode('utf-8', errors='replace').split('\n')
    if a.lines:
        lo, hi = (int(x) for x in a.lines.split(':'))
        lines = lines[lo - 1:hi]
    text = '\n'.join(lines)
    out = dict(source=a.source, name=a.name, url=a.url, local_path=a.path, read_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
               file_bytes=len(raw), file_sha256=hashlib.sha256(raw).hexdigest(), line_range=a.lines, git_origin=git(a.path, 'config', '--get', 'remote.origin.url'),
               git_commit=git(a.path, 'rev-parse', 'HEAD'), note=a.note, text=text)
    os.makedirs(os.path.join(HERE, '..', 'evidence', a.source), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, '..', 'evidence', a.source, a.name + '.json'), 'w'), indent=1)
    print(f'{a.source}/{a.name}: {len(text)} characters, commit {out["git_commit"]}')


if __name__ == '__main__':
    main()
