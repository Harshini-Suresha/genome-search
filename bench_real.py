"""Real-genome measurements. Run after the other bench scripts:  python3 bench_real.py
Downloads two public NCBI RefSeq records (via the Biopython test suite on GitHub), writes FASTA copies
for the workbench into ../app/data, and adds a 'real' section to results.json."""
import json, pathlib, random, time, urllib.request
from collections import Counter, defaultdict
from Bio import SeqIO
HERE = pathlib.Path(__file__).resolve().parent; APPDATA = HERE.parent / 'app' / 'data'; APPDATA.mkdir(parents=True, exist_ok=True)
(HERE / 'data').mkdir(exist_ok=True)
BASE = 'https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/'
def get(name):
    p = HERE / 'data' / name
    if not p.exists(): urllib.request.urlretrieve(BASE + name, p)
    return p
rec = SeqIO.read(get('NC_000932.gb'), 'genbank'); seq = str(rec.seq).upper()
pl = SeqIO.read(get('NC_005816.fna'), 'fasta'); plasmid = str(pl.seq).upper()
gc = lambda s: (s.count('G') + s.count('C')) / len(s)
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rc = lambda s: ''.join(comp.get(c, 'N') for c in reversed(s))
rng = random.Random(3)
def iid(like):  # random sequence with the same length and base composition
    c = Counter(like); return ''.join(rng.choices('ACGT', [c[b] for b in 'ACGT'], k=len(like)))
R = dict(source='NCBI RefSeq records distributed with the Biopython test suite (public domain sequence data)',
         datasets=[dict(accession='NC_000932.1', name='Arabidopsis thaliana chloroplast', length=len(seq), gc=gc(seq)),
                   dict(accession='NC_005816.1', name='Yersinia pestis plasmid pPCP1', length=len(plasmid), gc=gc(plasmid))])

# genes: coding sequences with real coordinates
cds, used = [], Counter()
for f in rec.features:
    if f.type != 'CDS' or 'translation' not in f.qualifiers: continue
    g = (f.qualifiers.get('gene') or [f.qualifiers.get('locus_tag', ['cds'])[0]])[0]; used[g] += 1
    name = g if used[g] == 1 else f"{g}_{used[g]}"; nt = str(f.extract(rec.seq)).upper()
    cds.append(dict(name=name, start=int(f.location.start), end=int(f.location.end), strand=int(f.location.strand or 1), nt=nt, aa=f.qualifiers['translation'][0]))
for p, s in ((APPDATA / 'NC_000932.fasta', f">NC_000932 Arabidopsis thaliana chloroplast\n{seq}\n"), (APPDATA / 'NC_005816.fasta', f">NC_005816 Yersinia pestis plasmid pPCP1\n{plasmid}\n"),
             (APPDATA / 'NC_000932_cds.fasta', ''.join(f">{c['name']}\n{c['nt']}\n" for c in cds))): p.write_text(s)
R['genes'] = [dict(name=c['name'], start=c['start'], end=c['end'], strand=c['strand'], length=len(c['nt']), gc=round(gc(c['nt']), 4)) for c in cds]
R['gene_hist'] = dict(width=300, counts=[sum(1 for c in cds if b * 300 <= len(c['nt']) < (b + 1) * 300) for b in range(24)])

# GC along the genome
W, S = 2000, 1000
R['gc_window'] = dict(window=W, step=S, pos=list(range(0, len(seq) - W + 1, S)), gc=[round(gc(seq[i:i+W]), 4) for i in range(0, len(seq) - W + 1, S)])

# k-mer structure: real against composition-matched random
rnd = iid(seq); rows = []
for k in (6, 8, 10, 12, 14):
    row = dict(k=k)
    for label, s in (('real', seq), ('random', rnd)):
        c = Counter(s[i:i+k] for i in range(len(s) - k + 1)); n = len(s) - k + 1
        row[label] = dict(distinct=len(c), repeated_pct=100 * sum(v for v in c.values() if v > 1) / n, max_bucket=max(c.values()))
    rows.append(row)
R['kmer_real'] = rows

# inverted repeat: find the long stretch whose reverse complement occurs elsewhere in the genome
K = 24; fwd = defaultdict(list)
for i in range(len(seq) - K + 1): fwd[seq[i:i+K]].append(i)
hit = [i for i in range(len(seq) - K + 1) if any(abs(j - i) > K for j in fwd.get(rc(seq[i:i+K]), ()))]
clusters, cur = [], [hit[0], hit[0]]
for i in hit[1:]:
    if i - cur[1] <= 60: cur[1] = i
    else: clusters.append(cur); cur = [i, i]
clusters.append(cur); a, b = max(clusters, key=lambda c: c[1] - c[0]); partner = fwd[rc(seq[a:a+K])][0]
q = seq[a + 500:a + 560]; Kk = 11; idx = defaultdict(list)
for i in range(len(seq) - Kk + 1): idx[seq[i:i+Kk]].append(i)
find = lambda s: [p for p in idx.get(s[:Kk], ()) if seq.startswith(s, p)]
R['inverted_repeat'] = dict(start=a, end=b + K, length=b + K - a, partner_start=partner, query_len=len(q), forward_hits=len(find(q)), reverse_complement_hits=len(find(rc(q))))

# index and scan cost on real against random sequence of the same size
rows = []
for n in (10_000, 30_000, 100_000, len(seq)):
    row = dict(n=n)
    for label, s in (('real', seq[:n]), ('random', iid(seq[:n]))):
        t = time.perf_counter(); ix = defaultdict(list)
        for i in range(len(s) - 11 + 1): ix[s[i:i+11]].append(i)
        build = (time.perf_counter() - t) * 1000; qs = [s[p:p+24] for p in rng.sample(range(len(s) - 24), 200)]
        t = time.perf_counter()
        for _ in range(5):
            for qq in qs: [p for p in ix.get(qq[:11], ()) if s.startswith(qq, p)]
        look = (time.perf_counter() - t) / (5 * len(qs)) * 1000; best = 1e9
        for _ in range(3):
            t = time.perf_counter(); [i for i in range(len(s) - 24 + 1) if s[i:i+24] == qs[0]]; best = min(best, (time.perf_counter() - t) * 1000)
        row[label] = dict(build_ms=build, lookup_ms=look, scan_ms=best, max_bucket=max(len(v) for v in ix.values()))
    rows.append(row)
R['timing_real'] = rows

# proteins: all-against-all local alignment of the first 220 residues of 24 long proteins
def score(a, b, m=2, x=-1, g=-2):
    prev = [0] * (len(b) + 1); best = 0
    for ai in a:
        cur = [0]; c0 = 0
        for j, bj in enumerate(b, 1):
            v = max(0, prev[j-1] + (m if ai == bj else x), prev[j] + g, c0 + g); cur.append(v); c0 = v
            if v > best: best = v
        prev = cur
    return best
P = sorted((c for c in cds if len(c['aa']) >= 100), key=lambda c: -len(c['aa']))[:24]; seqs = [c['aa'][:220] for c in P]
selfs = [score(s, s) for s in seqs]; M = [[0.0] * len(P) for _ in P]
for i in range(len(P)):
    for j in range(i, len(P)):
        v = 1.0 if i == j else score(seqs[i], seqs[j]) / min(selfs[i], selfs[j]); M[i][j] = M[j][i] = round(v, 3)
R['protein_matrix'] = dict(names=[c['name'] for c in P], matrix=M, residues=220)
R['plasmid'] = dict(length=len(plasmid), gc=gc(plasmid))
d = json.load(open(HERE / 'results.json')); d['real'] = R; json.dump(d, open(HERE / 'results.json', 'w'), indent=1)
print('real section added;', len(cds), 'CDS; IR', R['inverted_repeat'])
