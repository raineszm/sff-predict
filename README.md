#

## Tools

At the moment this repo manages the python environment using [pixi](pixi.sh).

The data download and extraction proccess uses

- `curl` to download
- `gzip/gunzip` to extract the json file
- `jq` to prefilter the json file to somethinmg more manageable
- [`just`](just.systems) as a task runner
- `nbdime` for merge notebook changes

## Data Sources

Data were collected from the following sources:

- Goodreads book data circa 2017: <https://cseweb.ucsd.edu/~jmcauley/datasets/goodreads.html>
- Nebula award winners: <https://nebulas.sfwa.org/awards-by-year/>
- Hugo award winners: <https://www.thehugoawards.org/hugo-history/>
- World Fantasy award winners: <https://worldfantasy.org/awards/winners/>
- A larger list of major scifi/fantasy awards is available here: <http://www.sfadb.com/Awards_Directory>

## Project Brainstorming

The basic idea is to try break data into years and train a ranking model on each year independently.

### Features

- Bibliographic data
- Publisher information
- Author information
- World events/trends for the award year
- Proxies for sales information?

### Target

- Either a binary class classification problem (won an award or not)
- or a composite award score (how many awards/nominations did a book get)

### Model

Maybe a model that computes a ranking score which allows to pairwise compare two books in a year to see which is more likely to win awards.

Something like

- [https://lightgbm.readthedocs.io/en/stable/pythonapi/lightgbm.LGBMRanker.html)(LightGBM)
- [LambdaRank](https://proceedings.neurips.cc/paper_files/paper/2006/file/af44c4c56f385c43f2529f9b1b018f6a-Paper.pdf)
- [RankNet](https://www.microsoft.com/en-us/research/wp-content/uploads/2005/08/icml_ranking.pdf)

### Loss function / Validation

Measure [Discounted Cumulative Gain](https://en.wikipedia.org/wiki/Discounted_cumulative_gain).
