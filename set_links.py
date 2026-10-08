#!/usr/bin/env python3
"""Fill in every link in one go.

    python3 set_links.py --github yourname/genomequery --streamlit https://yourapp.streamlit.app --site https://genomequery.vercel.app

Replaces the YOUR_* placeholders in the website, the Streamlit app, both notebooks and the READMEs.
Pass only the options you have; run it again later to add the rest."""
import argparse, pathlib, sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--github', help='GitHub user/repo, e.g. yourname/genomequery')
ap.add_argument('--streamlit', help='URL of the deployed Streamlit app')
ap.add_argument('--site', help='URL of the deployed website (Vercel)')
a = ap.parse_args()
repl = {}
if a.github: repl['YOUR_GITHUB_USER/YOUR_REPO'] = a.github.strip('/ ').removeprefix('https://github.com/')
if a.streamlit: repl['YOUR_STREAMLIT_URL'] = a.streamlit.rstrip('/')
if a.site: repl['YOUR_SITE_URL'] = a.site.rstrip('/')
if not repl: ap.print_help(); sys.exit(1)
root = pathlib.Path(__file__).resolve().parent; n = 0
for f in root.rglob('*'):
    if f.suffix in ('.html', '.md', '.ipynb', '.py', '.toml') and f.is_file() and f.name != 'set_links.py' and '__pycache__' not in f.parts:
        s = f.read_text(); o = s
        for k, v in repl.items(): s = s.replace(k, v)
        if s != o: f.write_text(s); n += 1; print('updated', f.relative_to(root))
print(f"{n} file(s) updated.")
