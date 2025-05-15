PARENT_JSON := justfile_dir() + "/data/goodreads_books.json"

fetch-dataset:
    [[ -f data/goodreads_books.json ]] && echo "Data file already downloaded. Delete it to redownload" && exit 1
    curl -sL https://mcauleylab.ucsd.edu/public_datasets/gdrive/goodreads/goodreads_books.json.gz -o data/goodreads_books.json.gz
    gunzip data/goodreads_books.json.gz

extract_scifi_fantasy:
    jq --slurp '[.[] | select(any(.popular_shelves[]?.name; test("fantasy|scifi|science fiction"; "i")))]' {{PARENT_JSON}} > {{justfile_dir()}}/data/scifi_fantasy.json
