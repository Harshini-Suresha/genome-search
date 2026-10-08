"""GenomeQuery Workbench: a working model of the system. Load sequences, index and search them,
align pairs, schedule and execute a dependency graph on real workers, query the relational store,
and benchmark everything on your own machine."""
import os, time, json, pathlib
import altair as alt, numpy as np, pandas as pd, streamlit as st
from engine import seqtools as sq, align as al, workflow as wf, bench as bn
from engine.store import Store, SCHEMA

CONTACT_EMAIL = "harshinisuresha7@gmail.com"
DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

def bundled_datasets():
    """Real example datasets shipped in app/data (see datasets.json)."""
    try:
        cat = json.loads((DATA_DIR / "datasets.json").read_text())
        return cat.get("datasets", [])
    except Exception:
        return []

alt.data_transformers.disable_max_rows()
GITHUB, SITE_URL = "YOUR_GITHUB_USER/YOUR_REPO", "YOUR_SITE_URL"   # filled in by set_links.py
BLUE, ORANGE, GREEN, SLATE, PLUM = "#0b63b5", "#e07a1f", "#1c8c6e", "#5b6b7c", "#8a4fb0"
KIND_COLORS = alt.Scale(domain=['gc', 'kmer', 'align', 'motif', 'report'], range=[BLUE, GREEN, ORANGE, PLUM, SLATE])
st.set_page_config(page_title="GenomeQuery Workbench", page_icon="🧬", layout="wide")
st.markdown("<style>.block-container{padding-top:2rem;max-width:1250px}h1,h2,h3{letter-spacing:-.01em}</style>", unsafe_allow_html=True)

def show(chart):
    try: st.altair_chart(chart, width="stretch")
    except TypeError: st.altair_chart(chart, use_container_width=True)
def tbl(df, **kw):
    try: st.dataframe(df, hide_index=True, width="stretch", **kw)
    except TypeError: st.dataframe(df, hide_index=True, use_container_width=True, **kw)

ss = st.session_state
def load(seqs, label):
    ss.store.reset(); ss.store.add_sequences(label, seqs, sq.kind_of, sq.gc)
    ss.seqs, ss.ver = seqs, ss.get('ver', 0) + 1
    for k in ('wf_res', 'idx_key', 'trie_key', 'sweep'): ss.pop(k, None)
if 'store' not in ss:
    ss.store, ss.seqs, ss.ver, ss.bench = Store(), {}, 0, {}
    load(sq.synthetic(), 'synthetic')
names = list(ss.seqs)

def get_index(k):
    if ss.get('idx_key') != (ss.ver, k): ss.idx = sq.build_index(ss.seqs, k); ss.idx_key = (ss.ver, k)
    return ss.idx

with st.sidebar:
    st.title("🧬 GenomeQuery")
    st.caption("Workbench: a working model of the system, not a mirror of the website. Everything here runs real code on your data and your machine.")
    st.metric("Sequences loaded", len(names)); st.metric("Residues", f"{sum(map(len, ss.seqs.values())):,}")
    st.metric("CPUs on this host", os.cpu_count())
    st.caption("Tabs follow the system: data, search, alignment, workflow engine, database, benchmarks.")
    st.caption(f"Contact: {CONTACT_EMAIL}")
    ok = lambda v: not v.startswith("YOUR_")
    if ok(SITE_URL): st.link_button("Project website", SITE_URL)
    if ok(GITHUB):
        st.link_button("Open in Colab", f"https://colab.research.google.com/github/{GITHUB}/blob/main/notebooks/genomequery_colab.ipynb")
        st.link_button("Source on GitHub", f"https://github.com/{GITHUB}")

tabs = st.tabs(["1 · Data", "2 · Search", "3 · Align", "4 · Workflow engine", "5 · Database", "6 · Benchmark"])

# ------------------------------------------------------------------ 1 DATA
with tabs[0]:
    left, right = st.columns([1, 1.3])
    with left:
        st.subheader("Load sequences")
        src = st.radio("Source", ["Generate synthetic", "Real example", "Upload FASTA", "Paste FASTA"], horizontal=True, key="src")
        if src == "Real example":
            ds = bundled_datasets()
            if not ds:
                st.info("No bundled datasets found in app/data. Use synthetic or upload instead.")
            else:
                labels = [f"{d['name']} ({d['length']:,} bp, GC {d['gc']})" for d in ds]
                pick = st.selectbox("Dataset", labels, key="greal")
                d = ds[labels.index(pick)]
                st.caption(f"{d['organism']} · {d['accession']} · {d['description']}")
                if st.button("Load dataset", key="grealbtn", type="primary"):
                    txt = (DATA_DIR / d['file']).read_text()
                    seqs = sq.parse_fasta(txt)
                    if seqs: load(seqs, d['file']); st.rerun()
                    else: st.error("Could not parse that dataset file.")
        elif src == "Generate synthetic":
            n = st.slider("Sequences", 4, 60, 12, key="gn"); L = st.slider("Length", 500, 6000, 2000, 250, key="gl")
            g = st.slider("GC content", .30, .70, .50, .01, key="gg"); rel = st.slider("Fraction related to sequence 1", 0.0, 1.0, .4, .1, key="grel")
            motif = st.text_input("Planted motif", "GATTACAGATTACA", key="gm")
            if st.button("Generate and load", key="gbtn", type="primary"): load(sq.synthetic(n, L, g, motif.upper(), rel), 'synthetic'); st.rerun()
        elif src == "Upload FASTA":
            f = st.file_uploader("FASTA file", type=["fa", "fasta", "fna", "txt"], key="gup")
            if f is not None and st.button("Load file", key="gfile", type="primary"):
                d = sq.parse_fasta(f.getvalue().decode(errors='ignore'))
                if d: load(d, f.name); st.rerun()
                else: st.error("No sequences found in that file.")
        else:
            txt = st.text_area("FASTA text", ">a\nGATTACAGATTACA\n>b\nGATTACAGATAACA", height=140, key="gtxt")
            if st.button("Load text", key="gpaste", type="primary"):
                d = sq.parse_fasta(txt)
                if d: load(d, 'pasted'); st.rerun()
                else: st.error("No sequences found.")
        st.caption("Every load is stored in the relational schema on the Database tab.")
    with right:
        df = ss.store.df("SELECT name, kind, length, ROUND(gc,3) AS gc FROM sequences ORDER BY seq_id")
        c1, c2 = st.columns(2)
        with c1: show(alt.Chart(df).mark_bar(color=BLUE).encode(x=alt.X('length:Q', bin=alt.Bin(maxbins=15), title='Length'), y=alt.Y('count()', title='Sequences')).properties(height=190, title='Length distribution'))
        with c2: show(alt.Chart(df).mark_circle(size=70, color=ORANGE).encode(x=alt.X('length:Q', title='Length'), y=alt.Y('gc:Q', scale=alt.Scale(zero=False), title='GC'), tooltip=['name', 'length', 'gc']).properties(height=190, title='GC against length'))
        tbl(df, height=240)
        ds = bundled_datasets()
        if ds:
            st.markdown("**Bundled real datasets (app/data)**")
            tbl(pd.DataFrame([dict(dataset=d['name'], accession=d['accession'], sequences=d['sequences'], length=d['length'], gc=d['gc']) for d in ds]))

# ------------------------------------------------------------------ 2 SEARCH
with tabs[1]:
    c1, c2, c3 = st.columns([1, 1.4, 1.6])
    k = c1.slider("k-mer length", 4, 16, 8, key="k")
    mode = c2.radio("Search type", ["Exact", "Mismatches (seed and vote)", "Prefix (trie)"], key="smode")
    qraw = c3.text_input("Query (sequence, or identifier prefix for the trie)", "GATTACAGATTACA", key="q"); q = qraw.upper().strip()
    idx, bsecs, bmb = get_index(k); stt = idx.stats()
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Distinct k-mers", f"{stt['distinct']:,}"); m2.metric("Index entries", f"{stt['entries']:,}")
    m3.metric("Build time", f"{bsecs*1000:.1f} ms"); m4.metric("Peak build memory", f"{bmb:.1f} MB")
    def hitmap(rows):
        d = pd.DataFrame(rows, columns=['sequence', 'position'] + (['mismatches'] if len(rows and rows[0]) > 2 else []) + ([] if len(rows and rows[0]) <= 3 else ['votes']))
        lens = {n: len(s) for n, s in ss.seqs.items()}
        base = alt.Chart(pd.DataFrame(dict(sequence=list(lens), end=list(lens.values())))).mark_rule(color="#c9d8e4", size=6).encode(y=alt.Y('sequence:N', sort=None, title=None), x=alt.value(0), x2='end:Q')
        pts = alt.Chart(d).mark_tick(color=ORANGE, thickness=3, size=16).encode(y=alt.Y('sequence:N', sort=None), x=alt.X('position:Q', title='Position'), tooltip=list(d.columns))
        return base + pts
    if len(q) < k and mode != "Prefix (trie)": st.warning(f"The query needs at least {k} bases (the k-mer length).")
    elif mode == "Exact":
        t = time.perf_counter(); hits = idx.exact(q); t_idx = (time.perf_counter() - t) * 1000
        t = time.perf_counter(); ref = [(n, p) for n, s in ss.seqs.items() for p in sq.naive_find(s, q)]; t_scan = (time.perf_counter() - t) * 1000
        a, b, c = st.columns(3); a.metric("Matches", len(hits)); b.metric("Index lookup", f"{t_idx:.3f} ms"); c.metric("Scan of every sequence", f"{t_scan:.3f} ms", f"{t_scan/max(t_idx,1e-6):.0f}× slower", delta_color="off")
        st.caption("Index and scan agree on every match." if sorted(hits) == sorted(ref) else "Index and scan disagree: this is a bug.")
        if hits:
            show(hitmap(hits).properties(height=max(120, 18 * len(names)), title="Where the query matches"))
            tbl(pd.DataFrame([(n, p, ss.seqs[n][max(0, p-8):p] + '[' + ss.seqs[n][p:p+len(q)] + ']' + ss.seqs[n][p+len(q):p+len(q)+8]) for n, p in hits], columns=['sequence', 'position', 'context']))
        else: st.info("No exact matches. Try a shorter query or the mismatch search.")
    elif mode == "Mismatches (seed and vote)":
        d = st.slider("Mismatches allowed", 0, 6, 2, key="dmm"); t = time.perf_counter(); rows = idx.approx(q, d); ms = (time.perf_counter() - t) * 1000
        a, b, c = st.columns(3); a.metric("Matches within distance", len(rows)); b.metric("Search time", f"{ms:.2f} ms"); c.metric("Finds every match?", "Guaranteed" if idx.guaranteed(len(q), d) else "Not guaranteed")
        st.caption(f"Each query k-mer votes for a start position and candidates are verified. A guarantee holds when {len(q)-k+1} k-mers exceed {d}×{k} = {d*k}, since {d} mismatches break at most {d*k} of them. Shorten k to widen the guarantee.")
        if rows:
            show(hitmap(rows).properties(height=max(120, 18 * len(names)), title="Where the query matches"))
            tbl(pd.DataFrame(rows, columns=['sequence', 'position', 'mismatches', 'seed votes']))
        else: st.info("No matches within that distance.")
    else:
        by = st.radio("Search by", ["Start of sequence", "Sequence identifier"], horizontal=True, key="tby")
        depth = st.slider("Characters indexed per sequence", 8, 60, 24, key="tdepth") if by == "Start of sequence" else 0
        key = (ss.ver, by, depth)
        if ss.get('trie_key') != key:
            tr = sq.Trie()
            for n, s in ss.seqs.items(): tr.insert(s[:depth] if depth else n, n)
            ss.trie, ss.trie_key = tr, key
        tr = ss.trie; term = q if by == "Start of sequence" else qraw.strip(); t = time.perf_counter(); res = tr.prefix(term); ms = (time.perf_counter() - t) * 1000
        a, b, c = st.columns(3); a.metric("Sequences with that prefix", len(res)); b.metric("Trie nodes", f"{tr.nodes:,}"); c.metric("Lookup", f"{ms:.3f} ms")
        per = []; level = [tr.root]
        for dpt in range(1, 1 + max(1, min(depth or 20, 60))):
            level = [v for nd in level for kk, v in nd.items() if kk != '$']; per.append((dpt, len(level)))
            if not level: break
        show(alt.Chart(pd.DataFrame(per, columns=['depth', 'nodes'])).mark_line(point=True, color=BLUE).encode(x=alt.X('depth:Q', title='Depth'), y=alt.Y('nodes:Q', title='Nodes at depth')).properties(height=200, title='How the trie branches'))
        if res: st.write(", ".join(sorted(res)[:60]) + (" ..." if len(res) > 60 else ""))
        else: st.info("No sequence starts with that prefix. Identifiers look like seq_001; sequence queries use A, C, G and T.")
    with st.expander("K-mer bucket sizes: how many positions share a k-mer"):
        bs = pd.DataFrame(sorted(idx.bucket_sizes().items()), columns=['positions per k-mer', 'distinct k-mers'])
        show(alt.Chart(bs).mark_bar(color=PLUM).encode(x=alt.X('positions per k-mer:O'), y=alt.Y('distinct k-mers:Q', scale=alt.Scale(type='symlog')), tooltip=list(bs.columns)).properties(height=220))
        st.caption("Short k makes big buckets, which means many candidates to verify; long k makes small buckets but is fragile to mutation.")

# ------------------------------------------------------------------ 3 ALIGN
with tabs[2]:
    a1, a2, a3 = st.columns(3)
    an = a1.selectbox("Sequence A", names, key="aa"); bn_ = a2.selectbox("Sequence B", names, index=min(1, len(names) - 1), key="ab")
    amode = a3.radio("Mode", ["local", "global"], horizontal=True, key="amode")
    A_, B_ = ss.seqs[an], ss.seqs[bn_]; maxL = max(10, min(500, len(A_), len(B_)))
    s1, s2, s3, s4, s5 = st.columns(5)
    match = s1.slider("Match", 1, 5, 2, key="am"); mism = s2.slider("Mismatch", -5, 0, -1, key="ax"); gap = s3.slider("Gap", -6, -1, -2, key="ag")
    Lr = int(s4.number_input("Region length", 10, maxL, min(200, maxL), key="aL")); oa = int(s5.number_input("Start in A", 0, max(0, len(A_) - Lr), 0, key="aoa"))
    ob = int(st.number_input("Start in B", 0, max(0, len(B_) - Lr), 0, key="aob"))
    t = time.perf_counter(); r = al.align(A_[oa:oa+Lr], B_[ob:ob+Lr], amode, match, mism, gap); ms = (time.perf_counter() - t) * 1000
    m1, m2, m3, m4, m5 = st.columns(5); m1.metric("Score", r['score']); m2.metric("Identity", f"{r['identity']*100:.0f}%"); m3.metric("Gaps", r['gaps']); m4.metric("Aligned length", r['length']); m5.metric("Time", f"{ms:.0f} ms", f"{(Lr+1)**2:,} cells", delta_color="off")
    blocks = [f"A {oa + r['start'][0] + i + 1:>5} {r['a'][i:i+70]}\n        {r['mid'][i:i+70]}\nB {ob + r['start'][1] + i + 1:>5} {r['b'][i:i+70]}" for i in range(0, max(1, len(r['mid'])), 70)]
    st.code("\n\n".join(blocks) if r['mid'] else "No positive-scoring local alignment with these scores.", language=None)
    H = r['H']; f = max(1, -(-max(H.shape) // 90))
    P = H if f == 1 else H[:(H.shape[0] // f) * f, :(H.shape[1] // f) * f].reshape(H.shape[0] // f, f, H.shape[1] // f, f).max(axis=(1, 3))
    ii, jj = np.indices(P.shape); heat = alt.Chart(pd.DataFrame(dict(i=ii.ravel(), j=jj.ravel(), score=P.ravel()))).mark_rect().encode(
        x=alt.X('j:O', axis=None), y=alt.Y('i:O', axis=None), color=alt.Color('score:Q', scale=alt.Scale(scheme='blues'), legend=alt.Legend(title='Score')), tooltip=['i', 'j', 'score'])
    pth = pd.DataFrame(dict(i=[min(p[0] // f, P.shape[0] - 1) for p in r['path']], j=[min(p[1] // f, P.shape[1] - 1) for p in r['path']]))
    cA, cB = st.columns(2)
    with cA:
        st.markdown("**Scoring matrix with traceback**" + (f" (pooled {f}×{f})" if f > 1 else ""))
        show((heat + alt.Chart(pth).mark_circle(color=ORANGE, size=28).encode(x='j:O', y='i:O')).properties(height=360))
    with cB:
        w = st.slider("Dot plot word size", 5, 14, 9, key="dw"); X, Y = ss.seqs[an][:3000], ss.seqs[bn_][:3000]; tb = {}
        for j in range(len(Y) - w + 1): tb.setdefault(Y[j:j+w], []).append(j)
        pts = [(i, j) for i in range(len(X) - w + 1) for j in tb.get(X[i:i+w], ())][:15000]
        st.markdown("**Dot plot: shared words along both sequences**")
        show(alt.Chart(pd.DataFrame(pts, columns=['A position', 'B position'])).mark_circle(size=6, color=BLUE, opacity=.7).encode(
            x=alt.X('B position:Q'), y=alt.Y('A position:Q', scale=alt.Scale(reverse=True))).properties(height=360))
        st.caption("A diagonal line means the sequences are related over that stretch; scattered dots are chance matches.")

# ------------------------------------------------------------------ 4 WORKFLOW
with tabs[3]:
    cpu = os.cpu_count() or 1
    w1, w2, w3, w4 = st.columns(4)
    workers = w1.slider("Workers", 1, 16, min(4, max(2, cpu)), key="wk"); wmode = w2.radio("Pool type", ["threads", "processes"], horizontal=True, key="wm")
    cap = w3.slider("Resource budget (units)", 1, 24, min(24, workers + 2), key="wc"); reps = w4.slider("Work per task", 10, 400, 100, 10, key="wr")
    x1, x2, x3, x4 = st.columns(4)
    wk_k = x1.slider("k for k-mer tasks", 4, 14, 8, key="wkk"); wmotif = x2.text_input("Motif for scan tasks", "GATTACA", key="wmo").upper()
    alen = x3.slider("Alignment length", 100, 600, 300, 50, key="wal"); prio = x4.checkbox("Longest alignments first", True, key="wpr")
    base = wf.build_pipeline(ss.seqs, wk_k, wmotif, reps, alen)
    if not prio:
        for t_ in base: t_.priority = 0
    edf = pd.DataFrame([dict(id=t_.id, kind=t_.kind, deps=", ".join(t_.deps), priority=t_.priority, cost=t_.cost) for t_ in base])
    with st.expander(f"Task graph: {len(base)} tasks (edit dependencies, priority or cost)"):
        st.caption("Change a dependency to point at a later task to create a cycle and watch the scheduler refuse to run it.")
        ed = st.data_editor(edf, num_rows="dynamic", hide_index=True, key=f"ed_{hash((ss.ver, wk_k, wmotif, reps, alen, prio))}",
                            column_config={'kind': st.column_config.SelectboxColumn(options=list(wf.RUNNERS), required=True), 'priority': st.column_config.NumberColumn(step=1), 'cost': st.column_config.NumberColumn(min_value=1, step=1)})
    first = next(iter(ss.seqs.values())); DEF = dict(gc=dict(seq=first, reps=reps), kmer=dict(seq=first, k=wk_k, reps=max(1, reps // 20)), align=dict(a=first[:alen], b=first[:alen]), motif=dict(seq=first, motif=wmotif, reps=reps), report={})
    bp = {t_.id: t_ for t_ in base}; tasks = []
    for _, row in ed.iterrows():
        if pd.isna(row['id']) or not str(row['id']).strip(): continue
        tid, kind = str(row['id']).strip(), row['kind'] if not pd.isna(row['kind']) else 'report'
        deps = [d.strip() for d in ('' if pd.isna(row['deps']) else str(row['deps'])).split(',') if d.strip()]
        tasks.append(wf.Task(tid, kind, bp[tid].payload if tid in bp and bp[tid].kind == kind else DEF[kind], deps, int(0 if pd.isna(row['priority']) else row['priority']), int(1 if pd.isna(row['cost']) else row['cost'])))
    try: order, level = wf.toposort(tasks); ok = True
    except wf.CycleError as e: ok = False; st.error(f"{e}. The scheduler cannot order a graph with a cycle.")
    except KeyError as e: ok = False; st.error(str(e))
    if ok:
        pos, cnt = {}, {}
        for tid in order: cnt[level[tid]] = cnt.get(level[tid], 0) + 1; pos[tid] = (level[tid], cnt[level[tid]])
        by = {t_.id: t_ for t_ in tasks}
        nd = pd.DataFrame([dict(task=t_.id, kind=t_.kind, x=pos[t_.id][0], y=pos[t_.id][1], priority=t_.priority, cost=t_.cost) for t_ in tasks])
        ed_ = pd.DataFrame([dict(x=pos[d][0], y=pos[d][1], x2=pos[t_.id][0], y2=pos[t_.id][1]) for t_ in tasks for d in t_.deps])
        g = alt.Chart(nd).mark_circle(size=90).encode(x=alt.X('x:Q', title='Dependency level', axis=alt.Axis(tickMinStep=1)), y=alt.Y('y:Q', axis=None), color=alt.Color('kind:N', scale=KIND_COLORS), tooltip=['task', 'kind', 'priority', 'cost'])
        if len(ed_): g = alt.Chart(ed_).mark_rule(color="#9bb4c8", opacity=.35).encode(x='x:Q', y='y:Q', x2='x2:Q', y2='y2:Q') + g
        show(g.properties(height=260, title=f"Dependency graph ({max(level.values())+1} levels)"))
        b1, b2 = st.columns([1, 1])
        if b1.button("Run workflow", type="primary", key="wrun"):
            with st.spinner("Running on real workers..."):
                res = wf.run_workflow(tasks, workers, wmode, cap); ss.wf_res = res
                rid = ss.store.record_run(res['runs'], {n: i for n, i in ss.store.seq_ids().items()}); ss.wf_run_id = rid
        if b2.button(f"Sweep 1, 2, 4, 8 workers ({wmode})", key="wsweep"):
            rows = []
            with st.spinner("Running the workflow at each worker count..."):
                for wv in (1, 2, 4, 8):
                    rr = wf.run_workflow(tasks, wv, wmode, max(cap, wv + 2)); rows.append(dict(workers=wv, makespan_ms=rr['makespan'], serial_ms=rr['serial']))
            ss.sweep = pd.DataFrame(rows)
        res = ss.get('wf_res')
        if res:
            m = st.columns(6); m[0].metric("Wall time", f"{res['makespan']:.0f} ms"); m[1].metric("Serial CPU work", f"{res['serial']:.0f} ms"); m[2].metric("Speedup", f"{res['speedup']:.2f}×")
            m[3].metric("Worker CPU busy", f"{res['utilization']*100:.0f}%"); m[4].metric("Workers used", res['lanes']); m[5].metric("Peak budget used", f"{res['peak_cost']} / {max(cap, 1)}")
            gd = pd.DataFrame([dict(task=r['id'], kind=r['kind'], worker=f"worker {r['lane']+1}", start=r['start'], end=r['end'], cpu_ms=r['cpu_ms'], critical=r['id'] in res['critical']) for r in res['runs']])
            bars = alt.Chart(gd).mark_bar(height=16).encode(x=alt.X('start:Q', title='Milliseconds since start'), x2='end:Q', y=alt.Y('worker:N', title=None), color=alt.Color('kind:N', scale=KIND_COLORS), tooltip=['task', 'kind', alt.Tooltip('start:Q', format='.1f'), alt.Tooltip('end:Q', format='.1f')])
            crit = alt.Chart(gd[gd.critical]).mark_bar(height=16, fill=None, stroke=ORANGE, strokeWidth=2).encode(x='start:Q', x2='end:Q', y='worker:N')
            show((bars + crit).properties(height=max(120, 34 * res['lanes']), title=f"Real execution timeline on {res['lanes']} {wmode[:-1]} worker(s). Orange outline: critical path"))
            st.caption(f"Critical path: {' → '.join(res['critical'][:6])}{' …' if len(res['critical']) > 6 else ''} ({res['critical_ms']:.0f} ms). This CPU-time chain is a floor: no worker count can finish faster than it. Saved to the database as run {ss.get('wf_run_id')}.")
            with st.expander("Task results"):
                tbl(pd.DataFrame([dict(task=r['id'], kind=r['kind'], worker=r['lane'] + 1, start_ms=round(r['start'], 1), end_ms=round(r['end'], 1), cpu_ms=round(r['cpu_ms'], 1), result=str({k: v for k, v in (r['result'] or {}).items()})) for r in sorted(res['runs'], key=lambda r: r['start'])]))
        sw = ss.get('sweep')
        if sw is not None:
            sw = sw.assign(ideal_ms=sw.makespan_ms.iloc[0] / sw.workers)
            long = sw.melt('workers', ['makespan_ms', 'ideal_ms'], 'series', 'ms').replace({'makespan_ms': 'Measured', 'ideal_ms': 'Perfect scaling'})
            show(alt.Chart(long).mark_line(point=True).encode(x=alt.X('workers:O', title='Workers'), y=alt.Y('ms:Q', title='Wall time (ms)'), color=alt.Color('series:N', scale=alt.Scale(range=[BLUE, SLATE]), title=None)).properties(height=240, title='Measured scaling on this machine'))
            st.caption(f"This host reports {cpu} CPU(s). With one CPU, extra workers cannot help; with more, processes usually scale where threads are limited by Python's global lock.")
    with st.expander("Synchronization lab: why shared state needs a lock"):
        t1, t2 = st.columns(2); nt = t1.slider("Threads", 2, 16, 8, key="sth"); inc = t2.slider("Increments per thread", 100, 3000, 1000, 100, key="sin")
        if st.button("Run both versions", key="sbtn"):
            a, e, ta = wf.counter_demo(nt, inc, False); b, _, tb_ = wf.counter_demo(nt, inc, True)
            c1, c2, c3 = st.columns(3); c1.metric("Expected total", f"{e:,}"); c2.metric("Without a lock", f"{a:,}", f"{a-e:,} lost updates", delta_color="inverse"); c3.metric("With a lock", f"{b:,}", f"{tb_/ta:.1f}× the time", delta_color="off")
            st.caption("Each thread reads the counter, yields, then writes back read+1. Without mutual exclusion, threads overwrite each other. The lock serialises the critical section, which costs time but never loses an update.")

# ------------------------------------------------------------------ 5 DATABASE
with tabs[4]:
    store = ss.store; d1, d2 = st.columns([1, 2])
    with d1:
        on = st.toggle("Indexes on", value=store.indexed, key="ixon")
        if on != store.indexed: store.set_indexes(on)
        counts = pd.DataFrame([dict(table=t_, rows=int(store.df(f"SELECT COUNT(*) n FROM {t_}").n[0])) for t_ in ('samples', 'sequences', 'jobs', 'executions', 'results')])
        tbl(counts)
        with st.expander("Schema"): st.code(SCHEMA.strip(), language="sql")
    with d2:
        PRE = {"Longest sequences": "SELECT name, length, ROUND(gc,3) AS gc FROM sequences ORDER BY length DESC LIMIT 10",
               "Sequences per sample": "SELECT s.name AS sample, COUNT(*) AS n, ROUND(AVG(q.length)) AS avg_length FROM sequences q JOIN samples s USING(sample_id) GROUP BY s.name",
               "Job timeline (latest run)": "SELECT j.task, e.worker, ROUND(e.started_ms,1) AS start_ms, ROUND(e.finished_ms-e.started_ms,1) AS duration_ms FROM jobs j JOIN executions e USING(job_id) WHERE j.run_id=(SELECT MAX(run_id) FROM jobs) ORDER BY e.started_ms LIMIT 50",
               "Total work per task kind": "SELECT substr(j.task,1,instr(j.task,':')-1) AS kind, COUNT(*) AS jobs, ROUND(SUM(e.finished_ms-e.started_ms)) AS total_ms FROM jobs j JOIN executions e USING(job_id) WHERE instr(j.task,':')>0 GROUP BY kind ORDER BY total_ms DESC",
               "Alignment scores by sequence": "SELECT q.name, r.score FROM results r JOIN jobs j USING(job_id) JOIN sequences q ON q.seq_id=r.seq_id WHERE r.detail='align' ORDER BY r.score DESC LIMIT 20"}
        ch = st.selectbox("Example query", list(PRE), key="dpre"); sql = st.text_area("SQL (read-only)", PRE[ch], height=110, key=f"sql_{ch}")
        if st.button("Run query", key="dbrun", type="primary"):
            s_ = sql.strip().rstrip(';')
            if not s_.lower().startswith(('select', 'with', 'explain')): st.error("Only SELECT, WITH and EXPLAIN statements are allowed here.")
            else:
                try:
                    out = store.df(s_); t_ms = store.timed(s_, 5); tbl(out); st.caption(f"{len(out)} row(s), {t_ms:.3f} ms (best of 5).")
                    st.markdown("**Query plan**"); st.code("\n".join(store.explain(s_)), language=None)
                except Exception as e: st.error(str(e))
    st.divider(); st.subheader("Index experiment at scale")
    n_rows = st.select_slider("Rows in a synthetic table", [10_000, 50_000, 100_000, 250_000, 500_000], 100_000, key="dbn")
    if st.button("Time queries with and without indexes", key="dbexp"):
        with st.spinner("Building the table and timing both plans..."): ss.exp = store.index_experiment(n_rows)
    if ss.get('exp') is not None:
        e = ss.exp; show(alt.Chart(e).mark_bar().encode(x=alt.X('indexed:N', title=None, sort=['No indexes', 'With indexes']), y=alt.Y('ms:Q', scale=alt.Scale(type='log'), title='Milliseconds per query (log)'), color=alt.Color('indexed:N', scale=alt.Scale(range=[SLATE, BLUE]), legend=None), column=alt.Column('query:N', title=None), tooltip=['query', 'indexed', alt.Tooltip('ms:Q', format='.3f')]).properties(width=240, height=240))
        sp = e.pivot(index='query', columns='indexed', values='ms'); st.caption(" · ".join(f"{q_}: {sp.loc[q_, 'No indexes']/sp.loc[q_, 'With indexes']:.0f}× faster with indexes" for q_ in sp.index))

# ------------------------------------------------------------------ 6 BENCHMARK
with tabs[5]:
    st.write("Run the benchmark suite here, on the machine hosting this app. Nothing is precomputed.")
    b1, b2, b3 = st.columns([2, 1, 1])
    sizes = b1.multiselect("Reference sizes (bases)", [10_000, 30_000, 100_000, 300_000, 1_000_000], [10_000, 30_000, 100_000, 300_000], key="bsz")
    qlen = b2.slider("Query length", 12, 100, 24, key="bq"); reps = b3.slider("Repeats (fastest kept)", 1, 7, 3, key="brep")
    if st.button("Run search benchmark", key="bsearch", type="primary") and sizes:
        with st.spinner("Timing scans and the index..."): ss.bench['search'] = bn.bench_search(sorted(sizes), qlen, 11, reps)
    if 'search' in ss.bench:
        d = ss.bench['search']; show(alt.Chart(d).mark_line(point=True).encode(x=alt.X('n:Q', scale=alt.Scale(type='log'), title='Reference length (bases)'), y=alt.Y('ms:Q', scale=alt.Scale(type='log'), title='Milliseconds'), color=alt.Color('method:N', title=None), tooltip=['method', 'n', alt.Tooltip('ms:Q', format='.4f')]).properties(height=320, title='Search cost against reference size'))
        pv = d.pivot(index='n', columns='method', values='ms'); tbl(pv.reset_index().round(4)); st.download_button("Download CSV", d.to_csv(index=False), "search_benchmark.csv", key="bdl1")
        st.caption("The C-implemented str.find scan is included on purpose: it shows how much of the index's advantage depends on comparing against Python-level loops. Each run also checks that all three methods return identical matches.")
    st.divider(); lens = st.multiselect("Alignment lengths", [50, 100, 200, 400, 800], [50, 100, 200, 400], key="bal")
    if st.button("Run alignment benchmark", key="balign") and lens:
        with st.spinner("Filling matrices..."): ss.bench['align'] = bn.bench_align(sorted(lens))
    if 'align' in ss.bench:
        d = ss.bench['align']; show(alt.Chart(d).mark_line(point=True).encode(x=alt.X('length:Q', title='Sequence length (both)'), y=alt.Y('ms:Q', title='Milliseconds'), color=alt.Color('method:N', title=None), tooltip=['method', 'length', alt.Tooltip('ms:Q', format='.1f')]).properties(height=300, title='Alignment time against length'))
        st.download_button("Download CSV", d.to_csv(index=False), "align_benchmark.csv", key="bdl2")

st.divider()
st.caption(f"GenomeQuery workbench · questions or data requests: {CONTACT_EMAIL}")
