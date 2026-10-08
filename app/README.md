# GenomeQuery Workbench (Streamlit)

A working model of the system, not a copy of the website. It runs real code on your data.

```
pip install -r requirements.txt
streamlit run app.py
```

| Tab | What it does |
|---|---|
| 1 Data | Generate synthetic sequences (with related families and a planted motif), upload FASTA, or paste FASTA. Everything is stored in the relational schema. |
| 2 Search | k-mer index; exact search checked against a full scan; mismatch search by seed-and-vote with a stated sensitivity guarantee; trie prefix search; hit map across sequences. |
| 3 Align | Smith-Waterman and Needleman-Wunsch with editable scores, traceback, pooled score heatmap and dot plot. |
| 4 Workflow engine | Editable task graph, Kahn topological sort with cycle detection, priority queue, resource budget, thread or process pool, real timeline, critical path, worker-count sweep, and a lock demonstration. Runs are saved to the database. |
| 5 Database | Read-only SQL console, query plans, index toggle, and an index experiment on up to 500,000 rows. |
| 6 Benchmark | Search and alignment benchmarks measured on the host, with CSV download. |

`engine/` holds the code (`seqtools.py`, `align.py`, `workflow.py`, `store.py`, `bench.py`).
Check it with `python3 tests/test_engine.py`.

## Notes
- The store uses SQLite so the app runs anywhere. The schema is written for PostgreSQL's model (keys, constraints, joins, transactions); moving it is a matter of swapping the connection in `engine/store.py`.
- Speedup is reported as CPU work divided by wall time, so on a single-core host it stays near or below 1x. Run on a multi-core machine to see processes scale.
- Process pools use `fork` where available, otherwise `spawn`.

## Hosting
Vercel serves the static site (`../index.html`) but cannot run a Streamlit server, which needs a long-lived process with websockets. Host the app on Streamlit Community Cloud, Hugging Face Spaces, Render or Fly.io, then run `python3 ../set_links.py --streamlit <url>` to put the link on the website.

## Bundled data
`data/` ships real public-domain sequences catalogued in `datasets.json`: the Arabidopsis chloroplast genome (NC_000932.1, 154 kb), its 85 coding sequences, a 10 kb demo subset, and the Yersinia plasmid pPCP1 (NC_005816.1, 9.6 kb). Pick them from “Real example” on the Data tab.

Contact: harshinisuresha7@gmail.com
