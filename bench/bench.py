"""GenomeQuery reference benchmarks. Run: python3 bench.py > results.json
Pure-Python reference implementations; swap in your own modules to re-measure."""
import random, time, json, sys, heapq, sqlite3, tracemalloc, platform, os, bisect
from collections import defaultdict, deque
random.seed(42)
def dna(n): return ''.join(random.choices('ACGT', k=n))
def timeit(f, rep=3):
    best = 1e9
    for _ in range(rep):
        t = time.perf_counter(); r = f(); best = min(best, time.perf_counter() - t)
    return best * 1000, r
def mem(f):
    tracemalloc.start(); r = f(); _, p = tracemalloc.get_traced_memory(); tracemalloc.stop(); return p / 1e6, r
R = {"env": {"python": platform.python_version(), "cpus": os.cpu_count(), "platform": platform.machine()}}

# 1. exact search: naive vs k-mer index
K = 11; M = 24
def naive(s, q):
    m = len(q); return [i for i in range(len(s) - m + 1) if s[i:i+m] == q]
def build(s, k=K):
    d = defaultdict(list)
    for i in range(len(s) - k + 1): d[s[i:i+k]].append(i)
    return d
def lookup(d, s, q):
    return [p for p in d.get(q[:K], ()) if s[p:p+len(q)] == q]
rows = []
for n in [10_000, 30_000, 100_000, 300_000, 1_000_000]:
    s = dna(n); p = random.randrange(n - M); q = s[p:p+M]
    tn, a = timeit(lambda: naive(s, q), 3)
    tb, d = timeit(lambda: build(s), 1)
    tl, b = timeit(lambda: lookup(d, s, q), 20)
    assert a == b
    mb, _ = mem(lambda: build(s))
    rows.append(dict(n=n, naive_ms=tn, build_ms=tb, lookup_ms=tl, index_mb=mb))
R["search"] = rows

# 1b. k sensitivity at n=300k
s = dna(300_000); ks = []
for k in [6, 8, 10, 12, 14]:
    tb, d = timeit(lambda: build(s, k), 1)
    q = s[1000:1024]
    t = time.perf_counter()
    for _ in range(200): d.get(q[:k], ())
    ks.append(dict(k=k, build_ms=tb, kmers=len(d), avg_bucket=(len(s) - k + 1) / len(d)))
R["ksweep"] = ks

# 2. prefix search: trie vs sorted list+bisect vs linear scan
class T:
    __slots__ = ("c", "e")
    def __init__(s): s.c = {}; s.e = 0
def tadd(root, w):
    for ch in w: root = root.c.setdefault(ch, T())
    root.e += 1
def tcount(root, p):
    for ch in p:
        root = root.c.get(ch)
        if root is None: return 0
    st = [root]; n = 0
    while st:
        x = st.pop(); n += x.e; st.extend(x.c.values())
    return n
rows = []
for N in [1_000, 5_000, 20_000, 50_000]:
    W = [dna(32) for _ in range(N)]; pre = W[N // 2][:10]
    root = T(); tb, _ = timeit(lambda: [tadd(root, w) for w in W], 1)
    SW = sorted(W)
    tt, a = timeit(lambda: tcount(root, pre), 20)
    ts, b = timeit(lambda: bisect.bisect_left(SW, pre + 'T' * 22) - bisect.bisect_left(SW, pre), 20)
    tl, c = timeit(lambda: sum(w.startswith(pre) for w in W), 5)
    mt, _ = mem(lambda: [tadd(T(), w) for w in W[:1]] and None)
    rows.append(dict(n=N, trie_us=tt*1000, bisect_us=ts*1000, scan_us=tl*1000, build_ms=tb))
R["prefix"] = rows

# 3. alignment: Smith-Waterman, full matrix vs two-row
def sw_full(a, b, m=2, x=-1, g=-2):
    H = [[0]*(len(b)+1) for _ in range(len(a)+1)]; best = 0
    for i in range(1, len(a)+1):
        ai = a[i-1]; Hi = H[i]; Hp = H[i-1]
        for j in range(1, len(b)+1):
            v = max(0, Hp[j-1] + (m if ai == b[j-1] else x), Hp[j] + g, Hi[j-1] + g)
            Hi[j] = v
            if v > best: best = v
    return best
def sw_two(a, b, m=2, x=-1, g=-2):
    prev = [0]*(len(b)+1); best = 0
    for i in range(1, len(a)+1):
        cur = [0]*(len(b)+1); ai = a[i-1]
        for j in range(1, len(b)+1):
            v = max(0, prev[j-1] + (m if ai == b[j-1] else x), prev[j] + g, cur[j-1] + g)
            cur[j] = v
            if v > best: best = v
        prev = cur
    return best
rows = []
for L in [100, 200, 400, 800, 1200]:
    a, b = dna(L), dna(L)
    tf, r1 = timeit(lambda: sw_full(a, b), 1); mf, _ = mem(lambda: sw_full(a, b))
    t2, r2 = timeit(lambda: sw_two(a, b), 1); m2, _ = mem(lambda: sw_two(a, b))
    assert r1 == r2
    rows.append(dict(len=L, full_ms=tf, two_ms=t2, full_mb=mf, two_mb=m2, cells=L*L))
R["align"] = rows

# 4. scheduler: Kahn + heap on random DAGs
def sched(n, edges_per):
    adj = defaultdict(list); indeg = [0]*n
    for v in range(1, n):
        for u in random.sample(range(max(0, v-50), v), min(edges_per, v, 50)): adj[u].append(v); indeg[v] += 1
    pri = [random.random() for _ in range(n)]
    h = [(pri[i], i) for i in range(n) if indeg[i] == 0]; heapq.heapify(h); order = []
    while h:
        _, u = heapq.heappop(h); order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0: heapq.heappush(h, (pri[v], v))
    assert len(order) == n
    return sum(len(x) for x in adj.values())
rows = []
for n in [1_000, 5_000, 20_000, 50_000, 100_000]:
    t, e = timeit(lambda: sched(n, 3), 3); rows.append(dict(n=n, edges=e, ms=t))
R["sched_scale"] = rows

# 4b. workflow: measured task durations -> list-scheduling makespan vs workers
def task_cost(kind):
    a, b = dna(kind), dna(kind); t = time.perf_counter(); sw_two(a, b); return time.perf_counter() - t
stages = [("fetch", 1, 150), ("qc", 6, 250), ("index", 6, 300), ("align", 12, 400), ("call", 6, 300), ("report", 1, 200)]
nodes = []; layer = []
for name, cnt, L in stages:
    cur = []
    for i in range(cnt):
        nodes.append(dict(id=len(nodes), name=f"{name}-{i}", stage=name, dur=task_cost(L), deps=[x for x in (layer[i % len(layer)::] if layer and cnt >= len(layer) else layer)][:2 if cnt > 1 else 99]))
        cur.append(len(nodes) - 1)
    layer = cur
def makespan(nodes, w):
    indeg = {n["id"]: len(n["deps"]) for n in nodes}; ch = defaultdict(list)
    for n in nodes:
        for d in n["deps"]: ch[d].append(n["id"])
    ready = [(-n["dur"], n["id"]) for n in nodes if indeg[n["id"]] == 0]; heapq.heapify(ready)
    run = []; t = 0; done = 0; busy = 0
    dur = {n["id"]: n["dur"] for n in nodes}
    while ready or run:
        while ready and len(run) < w:
            _, i = heapq.heappop(ready); heapq.heappush(run, (t + dur[i], i)); busy += dur[i]
        t, i = heapq.heappop(run)
        for c in ch[i]:
            indeg[c] -= 1
            if indeg[c] == 0: heapq.heappush(ready, (-dur[c], c))
    return t, busy
total = sum(n["dur"] for n in nodes)
R["workflow"] = dict(tasks=len(nodes), serial_ms=total*1000, nodes=[dict(id=n["id"], stage=n["stage"], dur_ms=n["dur"]*1000, deps=n["deps"]) for n in nodes],
    makespan=[dict(workers=w, ms=makespan(nodes, w)[0]*1000, util=makespan(nodes, w)[1]/(w*makespan(nodes, w)[0])) for w in [1, 2, 3, 4, 6, 8, 12]])

# 5. database (SQLite): index vs no index, join
rows = []
for N in [10_000, 100_000, 500_000]:
    con = sqlite3.connect(":memory:"); c = con.cursor()
    c.execute("CREATE TABLE samples(sample_id INTEGER PRIMARY KEY, name TEXT)")
    c.execute("CREATE TABLE seqs(seq_id INTEGER PRIMARY KEY, sample_id INTEGER REFERENCES samples, length INT, gc REAL)")
    c.executemany("INSERT INTO samples VALUES(?,?)", [(i, f"s{i}") for i in range(1000)])
    c.executemany("INSERT INTO seqs VALUES(?,?,?,?)", [(i, random.randrange(1000), random.randrange(200, 5000), random.random()) for i in range(N)])
    con.commit()
    qry = "SELECT COUNT(*), AVG(gc) FROM seqs WHERE sample_id=?"
    join = "SELECT s.name, COUNT(*) FROM seqs q JOIN samples s USING(sample_id) WHERE q.length>4900 GROUP BY s.name"
    def rep(sql, a=()):
        t = time.perf_counter()
        for i in range(20): c.execute(sql, a or (i*7 % 1000,) if "?" in sql else ()).fetchall()
        return (time.perf_counter() - t) / 20 * 1000
    t_no = rep(qry); j_no = rep(join)
    c.execute("CREATE INDEX i1 ON seqs(sample_id)"); c.execute("CREATE INDEX i2 ON seqs(length)")
    t_ix = rep(qry); j_ix = rep(join)
    rows.append(dict(rows=N, scan_ms=t_no, idx_ms=t_ix, join_scan_ms=j_no, join_idx_ms=j_ix))
R["db"] = rows
json.dump(R, sys.stdout, indent=1)
