"""More GenomeQuery measurements. Run after bench.py and bench_extra.py:  python3 bench_more.py
Adds to results.json: variance, qlen, kmem, dictscale, memcmp, heap, dbwrite."""
import json, random, time, sqlite3, heapq, bisect, os, tempfile, tracemalloc
from collections import defaultdict
random.seed(11)
R = json.load(open('results.json'))
def dna(n): return ''.join(random.choices('ACGT', k=n))
pc = time.perf_counter
def build(s, k):
    d = defaultdict(list)
    for i in range(len(s) - k + 1): d[s[i:i+k]].append(i)
    return d
def mem(f):
    tracemalloc.start(); r = f(); _, p = tracemalloc.get_traced_memory(); tracemalloc.stop(); return p / 1e6, r

# 1. run-to-run variation: 30 repeats, 1M-base reference, 24-base query
K, M = 11, 24; s = dna(1_000_000); p = random.randrange(len(s) - M); q = s[p:p+M]; idx = build(s, K)
scan = []; look = []
for _ in range(30):
    t = pc(); [i for i in range(len(s) - M + 1) if s[i:i+M] == q]; scan.append((pc() - t) * 1000)
    t = pc()
    for _ in range(2000): [x for x in idx.get(q[:K], ()) if s[x:x+M] == q]
    look.append((pc() - t) / 2000 * 1000)
R['variance'] = dict(n=len(s), runs=30, scan_ms=scan, lookup_ms=look)

# 2. query length (300k reference)
s = dna(300_000); idx = build(s, K); rows = []
for m in (12, 24, 50, 100, 200):
    p = random.randrange(len(s) - m); q = s[p:p+m]
    best = 1e9
    for _ in range(3):
        t = pc(); [i for i in range(len(s) - m + 1) if s[i:i+m] == q]; best = min(best, (pc() - t) * 1000)
    t = pc()
    for _ in range(2000): [x for x in idx.get(q[:K], ()) if s[x:x+m] == q]
    rows.append(dict(m=m, scan_ms=best, lookup_ms=(pc() - t) / 2000 * 1000))
R['qlen'] = rows

# 3. index memory against k (300k reference)
R['kmem'] = [dict(k=k, mb=mem(lambda: build(s, k))[0]) for k in (6, 8, 10, 12, 14)]

# 4. hash-table lookup cost against table size
rows = []
for n in (1_000, 10_000, 100_000, 1_000_000):
    d = {i * 7919: i for i in range(n)}; keys = [random.randrange(n) * 7919 for _ in range(200_000)]
    best = 1e9
    for _ in range(3):
        t = pc()
        for k in keys: d[k]
        best = min(best, (pc() - t) / len(keys) * 1e9)
    rows.append(dict(n=n, ns=best))
R['dictscale'] = rows

# 5. memory per structure: 20,000 sequences of length 32
W = [dna(32) for _ in range(20_000)]
def trie():
    r = {}
    for w in W:
        d = r
        for ch in w: d = d.setdefault(ch, {})
    return r
def kmers(): return build(''.join(W[:2000]), 11)
R['memcmp'] = [dict(name=n, mb=mem(f)[0]) for n, f in [
    ('List of strings', lambda: list(W)), ('Sorted list', lambda: sorted(W)),
    ('Hash table', lambda: {w: i for i, w in enumerate(W)}), ('Trie', trie)]]

# 6. priority queue: microseconds per push+pop pair at different heap sizes
rows = []
for n in (1_000, 10_000, 100_000, 1_000_000):
    h = [random.random() for _ in range(n)]; heapq.heapify(h); vals = [random.random() for _ in range(50_000)]
    t = pc()
    for v in vals: heapq.heappushpop(h, v) if False else (heapq.heappush(h, v), heapq.heappop(h))
    rows.append(dict(n=n, us=(pc() - t) / len(vals) * 1e6))
R['heap'] = rows

# 7. database write throughput (file-backed SQLite)
def run(n, indexes, per_row):
    path = tempfile.mktemp(suffix='.db'); con = sqlite3.connect(path); c = con.cursor()
    c.execute("CREATE TABLE seqs(seq_id INTEGER PRIMARY KEY, sample_id INT, length INT, gc REAL)")
    if indexes: c.execute("CREATE INDEX a ON seqs(sample_id)"); c.execute("CREATE INDEX b ON seqs(length)")
    con.commit(); data = [(i, random.randrange(1000), random.randrange(200, 5000), random.random()) for i in range(n)]
    t = pc()
    if per_row:
        for r in data: c.execute("INSERT INTO seqs VALUES(?,?,?,?)", r); con.commit()
    else:
        c.executemany("INSERT INTO seqs VALUES(?,?,?,?)", data); con.commit()
    dt = pc() - t; con.close(); os.remove(path); return n / dt
R['dbwrite'] = [dict(name='One commit per row', rows=2000, per_s=run(2000, False, True)),
                dict(name='One transaction', rows=100000, per_s=run(100000, False, False)),
                dict(name='One transaction, 2 indexes', rows=100000, per_s=run(100000, True, False))]
json.dump(R, open('results.json', 'w'), indent=1)
print('added', ['variance', 'qlen', 'kmem', 'dictscale', 'memcmp', 'heap', 'dbwrite'])
