"""Extra Genome Search measurements. Run after bench.py:  python3 bench_extra.py
Reads results.json and adds: kdist, trie_depth, align_dist, seed, plans."""
import json, random, sqlite3
from collections import Counter, defaultdict
random.seed(7)
R = json.load(open('results.json'))
def dna(n): return ''.join(random.choices('ACGT', k=n))
def mutate(a, r): return ''.join(random.choice([x for x in 'ACGT' if x != c]) if random.random() < r else c for c in a)
def sw(a, b, m=2, x=-1, g=-2):
    prev = [0]*(len(b)+1); best = 0
    for i in range(1, len(a)+1):
        cur = [0]*(len(b)+1); ai = a[i-1]
        for j in range(1, len(b)+1):
            v = max(0, prev[j-1] + (m if ai == b[j-1] else x), prev[j] + g, cur[j-1] + g)
            cur[j] = v
            if v > best: best = v
        prev = cur
    return best

# 1. how many positions does each k-mer have? (300k reference)
s = dna(300_000); out = []
for k in (8, 10):
    c = Counter(s[i:i+k] for i in range(len(s) - k + 1)); h = Counter(min(v, 12) for v in c.values())
    out.append(dict(k=k, hist=[h.get(i, 0) for i in range(1, 13)], mean=(len(s) - k + 1) / len(c), distinct=len(c)))
R['kdist'] = out

# 2. trie: nodes per depth, 20,000 random sequences of length 32
N, L = 20_000, 32; root = {}
for _ in range(N):
    d = root
    for ch in dna(L): d = d.setdefault(ch, {})
levels = [root]; counts = []
for _ in range(L):
    levels = [c for n in levels for c in n.values()]; counts.append(len(levels))
R['trie_depth'] = dict(n=N, length=L, nodes=counts)

# 3. local alignment scores: unrelated vs diverged copies (length 100, 150 pairs each)
classes = [('Unrelated', None), ('10% mutated', .10), ('30% mutated', .30), ('50% mutated', .50)]; res = []
for name, r in classes:
    sc = []
    for _ in range(150):
        a = dna(100); b = dna(100) if r is None else mutate(a, r); sc.append(sw(a, b))
    h = [0]*21
    for v in sc: h[min(v // 10, 20)] += 1
    res.append(dict(name=name, hist=h, mean=sum(sc)/len(sc)))
R['align_dist'] = res

# 4. seed sensitivity: does any k-mer of a mutated query hit the true locus? (100k reference, 50-base queries)
ref = dna(100_000); QL = 50; rates = [0, .02, .05, .10, .15, .20]; ks = [6, 8, 10, 12, 14, 16]
recov, cand = [], []
for k in ks:
    idx = defaultdict(list)
    for i in range(len(ref) - k + 1): idx[ref[i:i+k]].append(i)
    row = []
    for r in rates:
        ok = 0
        for _ in range(200):
            p = random.randrange(len(ref) - QL); q = mutate(ref[p:p+QL], r)
            if any((p + j) in set(idx.get(q[j:j+k], ())) for j in range(QL - k + 1)): ok += 1
        row.append(ok / 200 * 100)
    recov.append(row)
    cand.append(sum(len(idx.get(ref[p:p+k], ())) for p in random.sample(range(len(ref) - k), 300)) / 300)
R['seed'] = dict(ks=ks, rates=rates, recov=recov, cand=cand, query_len=QL, ref_len=len(ref), trials=200)

# 5. query plans before and after indexing
con = sqlite3.connect(':memory:'); c = con.cursor()
c.execute("CREATE TABLE samples(sample_id INTEGER PRIMARY KEY, name TEXT)")
c.execute("CREATE TABLE seqs(seq_id INTEGER PRIMARY KEY, sample_id INTEGER REFERENCES samples, length INT, gc REAL)")
c.executemany("INSERT INTO samples VALUES(?,?)", [(i, f"s{i}") for i in range(1000)])
c.executemany("INSERT INTO seqs VALUES(?,?,?,?)", [(i, random.randrange(1000), random.randrange(200, 5000), random.random()) for i in range(20000)])
Q = {'filter': "SELECT COUNT(*), AVG(gc) FROM seqs WHERE sample_id=7",
     'join': "SELECT s.name, COUNT(*) FROM seqs q JOIN samples s USING(sample_id) WHERE q.length>4900 GROUP BY s.name"}
plan = lambda q: [r[3] for r in c.execute("EXPLAIN QUERY PLAN " + q)]
before = {k: plan(q) for k, q in Q.items()}
c.execute("CREATE INDEX i1 ON seqs(sample_id)"); c.execute("CREATE INDEX i2 ON seqs(length)")
R['plans'] = {k: dict(scan=before[k], index=plan(q)) for k, q in Q.items()}
R['plans']['sql'] = Q
json.dump(R, open('results.json', 'w'), indent=1)
print('added', [k for k in ('kdist', 'trie_depth', 'align_dist', 'seed', 'plans')])
