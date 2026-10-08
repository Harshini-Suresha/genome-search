"""Dependency-aware workflow scheduler and concurrent execution engine.

Graph + Kahn topological sort + priority queue decide WHAT may run; a thread or process pool
decides WHERE. A resource budget caps how much heavy work is in flight at once."""
import heapq, itertools, multiprocessing as mp, os, threading, time
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from .align import score_only

class CycleError(Exception):
    def __init__(self, nodes): super().__init__("Dependency cycle among: " + ", ".join(sorted(nodes))); self.nodes = nodes

@dataclass
class Task:
    id: str; kind: str; payload: dict = field(default_factory=dict)
    deps: list = field(default_factory=list); priority: int = 0; cost: int = 1

# ---- task bodies: top-level functions so they can be pickled for process pools
def _gc(p):
    s, n = p['seq'], 0
    for _ in range(p.get('reps', 1)): n = sum(1 for c in s if c in 'GC')
    return dict(score=round(n / max(1, len(s)), 4), position=len(s))
def _kmer(p):
    s, k = p['seq'], p.get('k', 8)
    for _ in range(p.get('reps', 1)): c = Counter(s[i:i+k] for i in range(len(s) - k + 1))
    top = c.most_common(1)[0] if c else ('', 0)
    return dict(score=len(c), position=0, top=top[0], top_count=top[1])
def _align(p):
    return dict(score=score_only(p['a'], p['b']), position=0)
def _motif(p):
    s, m, hits = p['seq'], p['motif'], 0
    for _ in range(p.get('reps', 1)):
        hits, i = 0, s.find(m)
        while i != -1 and m: hits += 1; i = s.find(m, i + 1)
    return dict(hits=hits, position=s.find(m))
def _report(p): return dict(score=0, position=0)
RUNNERS = dict(gc=_gc, kmer=_kmer, align=_align, motif=_motif, report=_report)

def _run(kind, payload):
    t0, c0 = time.time(), time.thread_time(); res = RUNNERS[kind](payload)
    return dict(result=res, start=t0, end=time.time(), cpu=time.thread_time() - c0,
                worker=os.getpid() if mp.current_process().name != 'MainProcess' else threading.get_ident())

def toposort(tasks):
    """Kahn's algorithm with a priority queue. Raises CycleError listing the tasks that can never start."""
    ids = {t.id for t in tasks}
    for t in tasks:
        for d in t.deps:
            if d not in ids: raise KeyError(f"{t.id} depends on unknown task {d}")
    indeg = {t.id: len(set(t.deps)) for t in tasks}; kids = defaultdict(list); by = {t.id: t for t in tasks}
    for t in tasks:
        for d in set(t.deps): kids[d].append(t.id)
    cnt = itertools.count(); h = [(-t.priority, next(cnt), t.id) for t in tasks if indeg[t.id] == 0]; heapq.heapify(h)
    order, level = [], {}
    while h:
        _, _, u = heapq.heappop(h); order.append(u); level[u] = max((level[d] + 1 for d in by[u].deps), default=0)
        for v in kids[u]:
            indeg[v] -= 1
            if indeg[v] == 0: heapq.heappush(h, (-by[v].priority, next(cnt), v))
    if len(order) != len(tasks): raise CycleError(ids - set(order))
    return order, level

def critical_path(runs):
    """Longest chain of CPU time through the dependency graph: a lower bound on wall time at any worker count."""
    by = {r['id']: r for r in runs}; dur = {r['id']: r['cpu_ms'] for r in runs}; best, prev = {}, {}
    for r in sorted(runs, key=lambda r: r['end']):
        ds = [d for d in r['deps'] if d in best]
        p = max(ds, key=lambda d: best[d], default=None); best[r['id']] = dur[r['id']] + (best[p] if p else 0); prev[r['id']] = p
    if not best: return [], 0.0
    end = max(best, key=best.get); path, cur = [], end
    while cur: path.append(cur); cur = prev[cur]
    return path[::-1], best[end]

def run_workflow(tasks, workers=4, mode='threads', capacity=None):
    order, level = toposort(tasks)
    by = {t.id: t for t in tasks}; indeg = {t.id: len(set(t.deps)) for t in tasks}; kids = defaultdict(list)
    for t in tasks:
        for d in set(t.deps): kids[d].append(t.id)
    cap = max(capacity or workers, max((t.cost for t in tasks), default=1))
    if mode == 'processes':
        ctx = mp.get_context('fork' if 'fork' in mp.get_all_start_methods() else 'spawn')
        ex = ProcessPoolExecutor(max_workers=workers, mp_context=ctx)
    else: ex = ThreadPoolExecutor(max_workers=workers)
    cnt = itertools.count(); heap = []
    def push(t): heapq.heappush(heap, (-t.priority, next(cnt), t.id))
    for t in tasks:
        if indeg[t.id] == 0: push(t)
    running, used, lanes, runs, peak = {}, 0, {}, [], 0
    t0 = time.time()
    try:
        while heap or running:
            deferred = []
            while heap and len(running) < workers:
                item = heapq.heappop(heap); t = by[item[2]]
                if used + t.cost > cap and running: deferred.append(item); continue
                running[ex.submit(_run, t.kind, t.payload)] = t; used += t.cost; peak = max(peak, used)
            for d in deferred: heapq.heappush(heap, d)
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for f in done:
                t = running.pop(f); used -= t.cost; r = f.result()
                lane = lanes.setdefault(r['worker'], len(lanes))
                runs.append(dict(id=t.id, kind=t.kind, lane=lane, start=(r['start'] - t0) * 1000, end=(r['end'] - t0) * 1000,
                                 cpu_ms=r['cpu'] * 1000, result=r['result'], deps=list(t.deps), priority=t.priority, cost=t.cost))
                for c in kids[t.id]:
                    indeg[c] -= 1
                    if indeg[c] == 0: push(by[c])
    finally: ex.shutdown(wait=True)
    make = max((r['end'] for r in runs), default=0.0); serial = sum(r['cpu_ms'] for r in runs)
    path, plen = critical_path(runs)
    return dict(runs=runs, makespan=make, serial=serial, speedup=serial / make if make else 1.0,
                utilization=sum(r['cpu_ms'] for r in runs) / (make * workers) if make else 0.0,
                critical=path, critical_ms=plen, lanes=len(lanes), peak_cost=peak, levels=level, order=order)

def build_pipeline(seqs, k=8, motif='GATTACA', reps=100, align_len=300, ref=None):
    """Per sequence: load (GC) -> k-mer spectrum, alignment to the reference, motif scan; then one report."""
    names = list(seqs); ref = ref or names[0]; tasks, leaves = [], []
    for n in names:
        s = seqs[n]
        tasks.append(Task(f"load:{n}", 'gc', dict(seq=s, reps=reps), [], priority=1))
        tasks.append(Task(f"kmer:{n}", 'kmer', dict(seq=s, k=k, reps=max(1, reps // 20)), [f"load:{n}"]))
        tasks.append(Task(f"align:{n}", 'align', dict(a=seqs[ref][:align_len], b=s[:align_len]), [f"load:{n}"], priority=len(s) // 100, cost=2))
        tasks.append(Task(f"motif:{n}", 'motif', dict(seq=s, motif=motif, reps=reps), [f"load:{n}"]))
        leaves += [f"kmer:{n}", f"align:{n}", f"motif:{n}"]
    tasks.append(Task("report", 'report', {}, leaves, priority=0))
    return tasks

def counter_demo(threads=8, increments=1500, use_lock=False):
    """Read-modify-write on a shared counter. Without a lock, updates are lost; with one, none are."""
    state, lock = {'n': 0}, threading.Lock()
    def work():
        for _ in range(increments):
            if use_lock:
                with lock: v = state['n']; time.sleep(0); state['n'] = v + 1
            else: v = state['n']; time.sleep(0); state['n'] = v + 1
    ts = [threading.Thread(target=work) for _ in range(threads)]; t = time.perf_counter()
    for x in ts: x.start()
    for x in ts: x.join()
    return state['n'], threads * increments, time.perf_counter() - t
