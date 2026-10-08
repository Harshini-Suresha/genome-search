# Launch checklist

Three things get launched, each in the place that suits it. Do them in this order, then run `set_links.py` once.

## 1. Put the project on GitHub
```
cd genome-search
git init && git add . && git commit -m "Genome Search"
# create an empty repo on github.com, then:
git branch -M main
git remote add origin https://github.com/<you>/<repo>.git
git push -u origin main
```
The Colab link only works once the repo is public and the notebooks are on `main`.

## 2. Website on Vercel
Import the repo at vercel.com/new. Framework preset: **Other**. No build command, no output directory (the site is `index.html` at the repo root). Vercel gives you a URL such as `https://<project>.vercel.app`.

## 3. Streamlit workbench on Streamlit Community Cloud
Vercel cannot run Streamlit (it needs a long-lived server with websockets), so host the app separately and free:
1. Go to share.streamlit.io and choose **Create app** from your GitHub repo.
2. Branch `main`, **Main file path `app/app.py`**. Dependencies come from `app/requirements.txt`.
3. Deploy. You get `https://<name>.streamlit.app`.

Other hosts that work: Hugging Face Spaces (Streamlit template), Render, Fly.io.

## 4. Notebooks
Nothing to deploy. `notebooks/genome-search-colab.ipynb` is self-contained and opens in Colab from the link; `notebooks/genome-search.ipynb` runs locally from the `notebooks/` folder.

## 5. Put every link everywhere
```
python3 set_links.py --github <you>/<repo> --streamlit https://<name>.streamlit.app --site https://<project>.vercel.app
git add . && git commit -m "Add links" && git push
```
This fills the links section of the website, the Streamlit sidebar, the Colab badge and the READMEs. Vercel and Streamlit redeploy on push.

Contact: harshinisuresha7@gmail.com
