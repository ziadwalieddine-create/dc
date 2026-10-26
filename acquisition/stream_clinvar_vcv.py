"""Stream the ClinVar VCV XML full release into Parquet tables without storing the 6 GB archive.

The weekly release dated 2026-10-18/19 (the same release as the VCF and tab-delimited files in data/external/clinvar) is read over HTTPS,
its MD5 is computed on the fly and compared with the .md5 published beside it, and every VariationArchive is parsed as it arrives.
Tables written to data/external/clinvar/vcv/ (all values as strings or integers, nothing interpreted):

  variation        one row per VariationArchive: identifiers, type, name, aggregate germline / somatic / oncogenicity classification with review
                   status, date and submission count, allele ids of the record, canonical SPDI, GRCh38 top-level HGVS, genes, cross-references
  allele_location  one row per SimpleAllele and assembly: Chr, accession, start/stop, inner/outer start/stop (copy-number and structural
                   alleles), VCF position and alleles, variant length
  scv              one row per submission (ClinicalAssertion): SCV accession, submitter, classification type and value, date, review status,
                   assertion method, origin, affected status, number tested, observation method types, conditions, citations, comment
  rcv              one row per RCV: conditions and the condition-level classification

The archive itself is not kept (the disk allowance is smaller than the file); the manifest records its size, MD5 (source-published and
computed), the XML ReleaseDate and the SHA-256 of every output table. Re-running replaces the tables.

  python3 stream_clinvar_vcv.py [--limit-records N]
"""
import argparse
import datetime
import hashlib
import io
import json
import os
import re
import sys
import time
import urllib.request
import zlib

import pyarrow as pa
import pyarrow.parquet as pq
from lxml import etree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import acquire  # noqa: E402

URL = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarVCVRelease_2026-1018.xml.gz'
TARGET = 'clinvar_vcv_xml_weekly_2026_1018'
OUT = os.environ.get('VCV_OUT') or os.path.join(ROOT, 'data', 'external', 'clinvar', 'vcv')
S = pa.string()
SCHEMAS = {
    'variation': pa.schema([(c, S) for c in ('variation_id', 'accession', 'version', 'record_type', 'variation_type', 'name', 'record_status', 'date_created', 'date_last_updated',
                                              'most_recent_submission', 'number_of_submissions', 'number_of_submitters', 'allele_ids', 'canonical_spdi', 'hgvs_grch38', 'genes', 'xrefs',
                                              'germline_review_status', 'germline_description', 'germline_date_last_evaluated', 'germline_submission_count',
                                              'somatic_review_status', 'somatic_description', 'somatic_date_last_evaluated', 'oncogenicity_review_status',
                                              'oncogenicity_description', 'oncogenicity_date_last_evaluated', 'release_date')]),
    'allele_location': pa.schema([(c, S) for c in ('variation_id', 'allele_id', 'allele_type', 'assembly', 'chr', 'accession', 'start', 'stop', 'inner_start', 'inner_stop', 'outer_start',
                                                   'outer_stop', 'display_start', 'display_stop', 'variant_length', 'position_vcf', 'reference_allele_vcf', 'alternate_allele_vcf', 'for_display')]),
    'scv': pa.schema([(c, S) for c in ('scv', 'version', 'variation_id', 'assertion_id', 'org_id', 'submitter', 'organization_category', 'classification_type', 'classification',
                                       'date_last_evaluated', 'review_status', 'submission_date', 'date_created', 'date_updated', 'contributes_to_aggregate', 'assertion_method',
                                       'origin', 'affected_status', 'number_tested', 'method_types', 'conditions', 'citations', 'comment', 'comment_length', 'record_status')]),
    'rcv': pa.schema([(c, S) for c in ('rcv', 'version', 'variation_id', 'title', 'conditions', 'germline_review_status', 'germline_description', 'germline_date_last_evaluated',
                                       'germline_submission_count', 'somatic_description', 'oncogenicity_description')]),
}


SEGMENT = int(os.environ.get("VCV_SEGMENT", 250_000))      # records per segment; a segment is the unit of resumption


class Writer:
    def __init__(self, name, segment):
        self.name, self.sch, self.buf, self.n = name, SCHEMAS[name], [], 0
        os.makedirs(os.path.join(OUT, name), exist_ok=True)
        self.path = os.path.join(OUT, name, f'part-{segment:04d}.parquet')
        self.w = pq.ParquetWriter(self.path + '.tmp', self.sch, compression='zstd')

    def add(self, row):
        self.buf.append(row)
        if len(self.buf) >= 50_000:
            self.flush()

    def flush(self):
        if self.buf:
            self.w.write_table(pa.Table.from_pylist([{k: (None if v is None else str(v)) for k, v in r.items()} for r in self.buf], schema=self.sch))
            self.n += len(self.buf)
            self.buf = []

    def close(self):
        self.flush()
        self.w.close()
        if self.n == 0:             # nothing was written for this segment: no empty file
            os.remove(self.path + '.tmp')
            return
        os.replace(self.path + '.tmp', self.path)


CLOSE = b'</VariationArchive>'


class Stream(io.RawIOBase):
    """Compressed HTTP body -> decompressed bytes; counts and hashes the compressed bytes. With skip=N the first N VariationArchive
    records are dropped (counted by their closing tag, not parsed) and the XML declaration and root start tag are re-emitted before
    the first kept record, so the parser sees a well-formed document that starts at record N+1."""

    def __init__(self, url, limit_bytes=None, skip=0):
        req = urllib.request.Request(url, headers={'User-Agent': acquire.UA})
        self.r = urllib.request.urlopen(req, timeout=120)
        self.d = zlib.decompressobj(zlib.MAX_WBITS | 16)
        self.md5, self.n, self.buf, self.done, self.limit = hashlib.md5(), 0, b'', False, limit_bytes
        self.last = time.time()
        self.skip, self.skipped, self.header = skip, 0, None
        self.pending = b''
        self.head_seen = b''

    def readable(self):
        return True

    def header_release_date(self):
        m = re.search(rb'ReleaseDate="([^"]+)"', self.head_seen)
        return m.group(1).decode() if m else None

    def _skip(self, data):
        """Drop whole records until `skip` closing tags have been seen; keep what follows the last dropped one."""
        buf = self.pending + data
        if self.header is None:
            i = buf.find(b'<ClinVarVariationRelease')
            j = buf.find(b'>', i) if i >= 0 else -1
            if j < 0:
                self.pending = buf
                return b''
            self.header = buf[:j + 1]
        pos = 0
        while self.skipped < self.skip:
            k = buf.find(CLOSE, pos)
            if k < 0:
                break
            self.skipped += 1
            pos = k + len(CLOSE)
        if self.skipped < self.skip:
            self.pending = buf[max(pos, len(buf) - len(CLOSE)):]
            return b''
        self.pending = b''
        return self.header + b'\n' + buf[pos:]

    def readinto(self, b):
        while len(self.buf) < len(b) and not self.done:
            chunk = self.r.read(1 << 22)
            if not chunk or (self.limit and self.n >= self.limit):
                self.buf += self.d.flush()
                self.done = True
                break
            self.md5.update(chunk)
            self.n += len(chunk)
            data = self.d.decompress(chunk)
            if len(self.head_seen) < 4000:
                self.head_seen += data[:4000]
            if self.skip and self.skipped < self.skip:
                data = self._skip(data)
            self.buf += data
            if time.time() - self.last > 60:
                print(f'  {self.n / 1e9:.2f} GB read', flush=True)
                self.last = time.time()
        k = min(len(b), len(self.buf))
        b[:k] = self.buf[:k]
        self.buf = self.buf[k:]
        return k


def t(e, path, default=None):
    x = e.find(path) if e is not None else None
    return (x.text or '').strip() if x is not None and x.text is not None else default


def a(e, name, default=None):
    return e.get(name, default) if e is not None else default


def cls_block(rec, tag):
    """Aggregate classification element (GermlineClassification / SomaticClinicalImpact / OncogenicityClassification) of a record."""
    c = rec.find(f'Classifications/{tag}') if rec is not None else None
    if c is None:
        return None, None, None, None
    d = c.find('Description')
    return t(c, 'ReviewStatus'), (d.text or '').strip() if d is not None and d.text else None, a(d, 'DateLastEvaluated'), a(d, 'SubmissionCount')


def parse_archive(va, W, release_date):
    vid = va.get('VariationID')
    rec = va.find('ClassifiedRecord')
    if rec is None:
        rec = va.find('IncludedRecord')
    gv = cls_block(rec, 'GermlineClassification')
    sv = cls_block(rec, 'SomaticClinicalImpact')
    ov = cls_block(rec, 'OncogenicityClassification')
    # top-level alleles of the record (not the alleles echoed inside ClinicalAssertion)
    alleles = []
    if rec is not None:
        for child in rec:
            if child.tag == 'SimpleAllele':
                alleles.append(child)
            elif child.tag in ('Haplotype', 'Genotype'):
                alleles.extend(child.iter('SimpleAllele'))
    first = alleles[0] if alleles else None
    genes, xrefs = [], []
    for al in alleles:
        for g in al.findall('GeneList/Gene'):
            genes.append(f"{g.get('Symbol')}:{g.get('GeneID')}")
        for x in al.findall('XRefList/XRef'):
            xrefs.append(f"{x.get('DB')}:{x.get('ID')}")
    hg38 = None
    if first is not None:
        for h in first.findall('HGVSlist/HGVS'):
            if h.get('Assembly') == 'GRCh38' and 'top-level' in (h.get('Type') or ''):
                hg38 = t(h, 'NucleotideExpression/Expression')
                break
    W['variation'].add(dict(variation_id=vid, accession=va.get('Accession'), version=va.get('Version'), record_type=va.get('RecordType'), variation_type=va.get('VariationType'),
                            name=va.get('VariationName'), record_status=t(va, 'RecordStatus'), date_created=va.get('DateCreated'), date_last_updated=va.get('DateLastUpdated'),
                            most_recent_submission=va.get('MostRecentSubmission'), number_of_submissions=va.get('NumberOfSubmissions'), number_of_submitters=va.get('NumberOfSubmitters'),
                            allele_ids='|'.join(x.get('AlleleID') for x in alleles), canonical_spdi=t(first, 'CanonicalSPDI'), hgvs_grch38=hg38, genes='|'.join(genes), xrefs='|'.join(xrefs),
                            germline_review_status=gv[0], germline_description=gv[1], germline_date_last_evaluated=gv[2], germline_submission_count=gv[3],
                            somatic_review_status=sv[0], somatic_description=sv[1], somatic_date_last_evaluated=sv[2], oncogenicity_review_status=ov[0],
                            oncogenicity_description=ov[1], oncogenicity_date_last_evaluated=ov[2], release_date=release_date))
    for al in alleles:
        for loc in al.findall('Location/SequenceLocation'):
            W['allele_location'].add(dict(variation_id=vid, allele_id=al.get('AlleleID'), allele_type=t(al, 'VariantType'), assembly=loc.get('Assembly'), chr=loc.get('Chr'),
                                          accession=loc.get('Accession'), start=loc.get('start'), stop=loc.get('stop'), inner_start=loc.get('innerStart'), inner_stop=loc.get('innerStop'),
                                          outer_start=loc.get('outerStart'), outer_stop=loc.get('outerStop'), display_start=loc.get('display_start'), display_stop=loc.get('display_stop'),
                                          variant_length=loc.get('variantLength'), position_vcf=loc.get('positionVCF'), reference_allele_vcf=loc.get('referenceAlleleVCF'),
                                          alternate_allele_vcf=loc.get('alternateAlleleVCF'), for_display=loc.get('forDisplay')))
    if rec is None:
        return
    for rcv in rec.findall('RCVList/RCVAccession'):
        conds = '|'.join(f"{c.get('DB')}:{c.get('ID')}:{(c.text or '').strip()}" for c in rcv.findall('ClassifiedConditionList/ClassifiedCondition'))
        rg, rs, ro = (cls_block_in(rcv.find('RCVClassifications'), k) for k in ('GermlineClassification', 'SomaticClinicalImpact', 'OncogenicityClassification'))
        W['rcv'].add(dict(rcv=rcv.get('Accession'), version=rcv.get('Version'), variation_id=vid, title=rcv.get('Title'), conditions=conds, germline_review_status=rg[0],
                          germline_description=rg[1], germline_date_last_evaluated=rg[2], germline_submission_count=rg[3], somatic_description=rs[1], oncogenicity_description=ro[1]))
    for ca in rec.findall('ClinicalAssertionList/ClinicalAssertion'):
        acc = ca.find('ClinVarAccession')
        cl = ca.find('Classification')
        ctype, cval = None, None
        if cl is not None:
            for k in ('GermlineClassification', 'SomaticClinicalImpact', 'OncogenicityClassification'):
                x = cl.find(k)
                if x is not None:
                    ctype, cval = k, (x.text or '').strip()
                    break
        oi = ca.find('ObservedInList/ObservedIn')
        samp = oi.find('Sample') if oi is not None else None
        methods = '|'.join(sorted({(m.text or '').strip() for m in ca.findall('ObservedInList/ObservedIn/Method/MethodType')}))
        origins = '|'.join(sorted({(m.text or '').strip() for m in ca.findall('ObservedInList/ObservedIn/Sample/Origin')}))
        affected = '|'.join(sorted({(m.text or '').strip() for m in ca.findall('ObservedInList/ObservedIn/Sample/AffectedStatus')}))
        ntested = '|'.join(sorted({(m.text or '').strip() for m in ca.findall('ObservedInList/ObservedIn/Sample/NumberTested')}))
        method = '|'.join(sorted({(x.text or '').strip() for x in ca.findall("AttributeSet/Attribute[@Type='AssertionMethod']")}))
        conds = '|'.join(sorted({f"{x.get('DB')}:{x.get('ID')}" for x in ca.findall('TraitSet/Trait/XRef')}))
        names = '|'.join(sorted({(x.text or '').strip() for x in ca.findall("TraitSet/Trait/Name/ElementValue[@Type='Preferred']")}))
        cites = '|'.join(sorted({f"{c.find('ID').get('Source')}:{(c.find('ID').text or '').strip()}" for c in ca.iter('Citation') if c.find('ID') is not None}))
        comment = t(cl, 'Comment') if cl is not None else None
        W['scv'].add(dict(scv=a(acc, 'Accession'), version=a(acc, 'Version'), variation_id=vid, assertion_id=ca.get('ID'), org_id=a(acc, 'OrgID'), submitter=a(acc, 'SubmitterName'),
                          organization_category=a(acc, 'OrganizationCategory'), classification_type=ctype, classification=cval, date_last_evaluated=a(cl, 'DateLastEvaluated'),
                          review_status=t(cl, 'ReviewStatus') if cl is not None else None, submission_date=ca.get('SubmissionDate'), date_created=ca.get('DateCreated'),
                          date_updated=ca.get('DateLastUpdated'), contributes_to_aggregate=ca.get('ContributesToAggregateClassification'), assertion_method=method, origin=origins,
                          affected_status=affected, number_tested=ntested, method_types=methods, conditions=(conds + ('#' + names if names else '')) or None, citations=cites or None,
                          comment=(comment[:2000] if comment else None), comment_length=len(comment) if comment else 0, record_status=t(ca, 'RecordStatus')))


def cls_block_in(parent, tag):
    c = parent.find(tag) if parent is not None else None
    if c is None:
        return None, None, None, None
    d = c.find('Description')
    return t(c, 'ReviewStatus'), (d.text or '').strip() if d is not None and d.text else None, a(d, 'DateLastEvaluated'), a(d, 'SubmissionCount')


PROGRESS = os.path.join(OUT, 'progress.json')


def load_progress():
    if os.path.exists(PROGRESS):
        d = json.load(open(PROGRESS))
        if d.get('url') == URL and d.get('segment') == SEGMENT:
            return d
    return dict(url=URL, segment=SEGMENT, done=[])


def save_progress(d):
    os.makedirs(OUT, exist_ok=True)
    with open(PROGRESS + '.tmp', 'w') as f:
        json.dump(d, f, indent=1)
    os.replace(PROGRESS + '.tmp', PROGRESS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit-records', type=int)
    ap.add_argument('--restart', action='store_true', help='discard finished segments')
    args = ap.parse_args()
    src_md5 = urllib.request.urlopen(urllib.request.Request(URL + '.md5', headers={'User-Agent': acquire.UA}), timeout=60).read().decode().split()[0]
    prog = dict(url=URL, segment=SEGMENT, done=[]) if args.restart else load_progress()
    if args.restart and os.path.isdir(OUT):
        import shutil
        shutil.rmtree(OUT)
    done = len(prog['done'])
    skip = done * SEGMENT
    st = Stream(URL, skip=skip)
    reader = io.BufferedReader(st, buffer_size=1 << 22)
    release_date, n, seg = prog.get('release_date'), skip, done
    if skip:
        print(f'resuming after {done} finished segments ({skip:,} records); the archive is read again from its start and the finished records are skipped, not parsed', flush=True)
    W = {k: Writer(k, seg) for k in SCHEMAS}
    t0 = time.time()
    head = None
    for ev, el in etree.iterparse(reader, events=('end',), tag='VariationArchive', huge_tree=True):
        if release_date is None:
            release_date = st.header_release_date() if hasattr(st, 'header_release_date') else None
        parse_archive(el, W, release_date)
        n += 1
        el.clear()
        while el.getprevious() is not None:
            del el.getparent()[0]
        if (n - skip) % SEGMENT == 0:
            for w in W.values():
                w.close()
            prog['done'].append(dict(segment=seg, records=SEGMENT, rows={k: w.n for k, w in W.items()}))
            prog['release_date'] = release_date
            save_progress(prog)
            print(f'segment {seg} written; {n:,} records, {st.n / 1e9:.2f} GB read, {time.time() - t0:.0f}s', flush=True)
            seg += 1
            W = {k: Writer(k, seg) for k in SCHEMAS}
        if args.limit_records and n - skip >= args.limit_records:
            break
    complete = not args.limit_records
    if complete:
        while st.read(1 << 22):
            pass
    for w in W.values():
        w.close()
    tail_rows = {k: w.n for k, w in W.items()}
    if tail_rows['variation'] and complete:
        prog['done'].append(dict(segment=seg, records=tail_rows['variation'], rows=tail_rows))
    prog['release_date'] = release_date
    save_progress(prog)
    got = st.md5.hexdigest()
    tables = {}
    for k in SCHEMAS:
        files = sorted(os.listdir(os.path.join(OUT, k)))
        tables[k] = dict(files=len(files), rows=sum(pq.ParquetFile(os.path.join(OUT, k, f)).metadata.num_rows for f in files if f.endswith('.parquet')),
                         bytes=sum(os.path.getsize(os.path.join(OUT, k, f)) for f in files),
                         sha256={f: hashlib.sha256(open(os.path.join(OUT, k, f), 'rb').read()).hexdigest() for f in files if f.endswith('.parquet')})
    ok = complete and got == src_md5
    rec = dict(id=TARGET, source='clinvar', url=URL, dest='clinvar/vcv', tier='A', kind='stream_derive', ok=ok, fetched_utc=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
               bytes=st.n, md5=got, source_md5=src_md5, md5_verified_against_source=got == src_md5, size_verified_against_source=False, release_date=release_date, records=n,
               tables={k: dict(files=v['files'], rows=v['rows'], bytes=v['bytes']) for k, v in tables.items()},
               sha256=hashlib.sha256(json.dumps({k: v['sha256'] for k, v in tables.items()}, sort_keys=True).encode()).hexdigest(),
               note='archive streamed and parsed, not stored; MD5 of the compressed stream compared with the published .md5')
    if not args.limit_records:
        acquire.update_manifest(rec)
        with open(os.path.join(HERE, 'manifest_files', TARGET + '.json'), 'w') as f:
            json.dump(dict(tables={k: v['sha256'] for k, v in tables.items()}, segments=prog['done']), f, indent=1)
    print(json.dumps({k: rec[k] for k in ('ok', 'bytes', 'md5', 'source_md5', 'release_date', 'records')}), json.dumps({k: v['rows'] for k, v in tables.items()}))


if __name__ == '__main__':
    main()
