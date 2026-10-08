# Genome Search

Integrated project: static site (`index.html`), Streamlit workbench (`app/` with `engine/`), benchmarks (`bench/` + `bench_real.py`), notebooks (`notebooks/`), and real example data (`app/data/` with `datasets.json`).

Contact: harshinisuresha7@gmail.com

Static site, no build step.

- `index.html` is the whole website (data is inlined in the `<script id="data">` block).
- `bench/` holds the benchmark programs that produced the data:
  `bench.py`, then `bench_extra.py`, then `bench_more.py` (run in that order from inside `bench/`).

Deploy on Vercel: import this folder with no framework and no build command, or run `vercel` inside it.

To refresh the numbers: run the three scripts, then paste the contents of `bench/results.json`
into the `<script id="data" type="application/json">` block of `index.html`.

## Workbench
`app/` is a Streamlit app that runs the system itself. `notebooks/` holds a Jupyter notebook (`genome-search.ipynb`) and a self-contained Colab notebook (`genome-search-colab.ipynb`). Launch steps are in `DEPLOY.md`; `set_links.py` puts every link into the site, app and notebooks.
