"""Run:  python3 tests/test_engine.py   (from the app folder)  or  pytest tests"""
import sys, random, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from engine import seqtools as sq, align as al, workflow as wf, bench as bn
from engine.store import Store

def test_search():
    seqs = sq.synthetic(8, 1500, seed=3); idx, _, _ = sq.build_index(seqs, 8)
    q = seqs['seq_001'][100:130]
    assert sorted(idx.exact(q)) == sorted((n, p) for n, s in seqs.items() for p in sq.naive_find(s, q))
    mq = sq.mutate(q, .06, random.Random(1), indel=0); d = sq.hamming(q, mq)
    assert any(r[0] == 'seq_001' and r[1] == 100 for r in idx.approx(mq, d))
    tr = sq.Trie()
    for n, s in seqs.items(): tr.insert(s[:12], n)
    assert 'seq_001' in tr.prefix(seqs['seq_001'][:6])
    assert sq.parse_fasta(">a x\nACGT\nAC\n>b\nGG") == {'a': 'ACGTAC', 'b': 'GG'}

def test_align():
    assert al.align('GGTTGACTA', 'TGTTACGG', 'local', 3, -3, -2)['score'] == 13
    assert al.score_only('GGTTGACTA', 'TGTTACGG', 3, -3, -2) == 13
    assert al.align('ACGT', 'ACGT', 'global')['identity'] == 1.0

def test_workflow():
    seqs = sq.synthetic(6, 1200, seed=3); T = wf.build_pipeline(seqs, reps=10, align_len=100)
    order, _ = wf.toposort(T); pos = {t: i for i, t in enumerate(order)}
    assert all(pos[d] < pos[t.id] for t in T for d in t.deps)
    try: wf.toposort([wf.Task('a', 'gc', deps=['b']), wf.Task('b', 'gc', deps=['a'])]); assert False
    except wf.CycleError: pass
    for mode in ('threads', 'processes'):
        res = wf.run_workflow(T, workers=3, mode=mode, capacity=4)
        end = {r['id']: r['end'] for r in res['runs']}; st = {r['id']: r['start'] for r in res['runs']}
        assert len(res['runs']) == len(T) and all(st[t.id] >= end[d] - 1 for t in T for d in t.deps)
        assert res['peak_cost'] <= 4
    n, exp, _ = wf.counter_demo(6, 300, False); n2, _, _ = wf.counter_demo(6, 300, True)
    assert n2 == exp and n < exp

def test_store_and_bench():
    seqs = sq.synthetic(5, 800, seed=3); st = Store(); st.add_sequences('S', seqs, sq.kind_of, sq.gc)
    res = wf.run_workflow(wf.build_pipeline(seqs, reps=5, align_len=80), workers=2)
    st.record_run(res['runs'], st.seq_ids()); assert st.df('SELECT COUNT(*) n FROM jobs').n[0] == len(res['runs'])
    assert 'SCAN' in st.explain("SELECT * FROM sequences WHERE sample_id=1")[0]
    st.set_indexes(True); assert 'SEARCH' in st.explain("SELECT * FROM sequences WHERE sample_id=1")[0]
    assert len(bn.bench_search([5000], reps=1)) == 4 and len(bn.bench_align([40])) == 2

if __name__ == '__main__':
    for f in (test_search, test_align, test_workflow, test_store_and_bench): f(); print('ok', f.__name__)
