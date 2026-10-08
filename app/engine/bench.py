"""Live benchmarks that run on the machine hosting the app."""
import random, time, tracemalloc
import pandas as pd
from .seqtools import KmerIndex, naive_find
from .align import align, score_only

def _best(f, reps):
    b = 1e9
    for _ in range(reps):
        t = time.perf_counter(); f(); b = min(b, time.perf_counter() - t)
    return b * 1000

def bench_search(sizes, qlen=24, k=11, reps=3, seed=1):
    rng = random.Random(seed); rows = []
    for n in sizes:
        s = ''.join(rng.choices('ACGT', k=n)); p = rng.randrange(n - qlen); q = s[p:p+qlen]
        def loop(): return [i for i in range(n - qlen + 1) if s[i:i+qlen] == q]
        tracemalloc.start(); t = time.perf_counter(); idx = KmerIndex(k); idx.add('s', s); build = (time.perf_counter() - t) * 1000
        peak = tracemalloc.get_traced_memory()[1] / 1e6; tracemalloc.stop()
        assert [x for _, x in idx.exact(q)] == loop() == naive_find(s, q)
        t = time.perf_counter()
        for _ in range(500): idx.exact(q)
        look = (time.perf_counter() - t) / 500 * 1000
        rows += [dict(n=n, method='Python loop scan', ms=_best(loop, reps), mb=None),
                 dict(n=n, method='str.find scan (C)', ms=_best(lambda: naive_find(s, q), reps), mb=None),
                 dict(n=n, method='K-mer index build (once)', ms=build, mb=peak),
                 dict(n=n, method='K-mer index lookup', ms=look, mb=None)]
    return pd.DataFrame(rows)

def bench_align(lengths, reps=1, seed=2):
    rng = random.Random(seed); rows = []
    for L in lengths:
        a, b = ''.join(rng.choices('ACGT', k=L)), ''.join(rng.choices('ACGT', k=L))
        rows.append(dict(length=L, method='Two-row score', ms=_best(lambda: score_only(a, b), reps)))
        if L <= 800: rows.append(dict(length=L, method='Full matrix + traceback', ms=_best(lambda: align(a, b), reps)))
    return pd.DataFrame(rows)
