import json

SCIFI_FANTASY_SHELVES = ["fantasy", "scifi", "science fiction"]
RATINGS_THRESHOLD = 1000
TAG_THRESHOLD = RATINGS_THRESHOLD * 0.05


def shelf_matches(shelf):
    return (
        shelf["name"] in SCIFI_FANTASY_SHELVES and int(shelf["count"]) >= TAG_THRESHOLD
    )


def is_scifi_fantasy(book):
    if not book["ratings_count"] or int(book["ratings_count"]) < RATINGS_THRESHOLD:
        return False
    return any(shelf_matches(shelf) for shelf in book["popular_shelves"])


with open("data/goodreads_books.json", "r") as f:
    with open("data/scifi_fantasy.json", "w") as out:
        for line in f:
            book = json.loads(line)
            if is_scifi_fantasy(book):
                out.write(json.dumps(book) + "\n")
