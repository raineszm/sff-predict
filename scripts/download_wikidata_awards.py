import sys
import httpx
import json

endpoint_url = "https://query.wikidata.org/sparql"

# This is sort of impenetrible without some reference to what
# the various entities and relations mean
# I reccomend pasting it into the online editor
# https://query.wikidata.org/
# which gives hover information for each entity and relation
query = """
SELECT DISTINCT ?title (GROUP_CONCAT(DISTINCT ?authorLabel; SEPARATOR=', ') as ?authors) ?year ?awardLabel
(GROUP_CONCAT(DISTINCT ?genreLabel; SEPARATOR=';') as ?genres)
(GROUP_CONCAT(DISTINCT ?work_id; SEPARATOR=';') as ?work_ids)
(GROUP_CONCAT(DISTINCT ?ol_id; SEPARATOR=';') as ?ol_ids)
(GROUP_CONCAT(DISTINCT ?isfdb_id; SEPARATOR=';') as ?isdfb_ids)
?status
WHERE {

  
  
  ## --- nominations -----------------------------------------------------------
  {
    ?item p:P1411 ?stmt .
    ?stmt ps:P1411 ?award ;
          pq:P585  ?date .
    BIND("nominated" AS ?status)
  }
  UNION
  ## --- wins ------------------------------------------------------------------
  {
    ?item p:P166  ?stmt .
    ?stmt ps:P166  ?award ;
          pq:P585  ?date .
    BIND("winner" AS ?status)
  }

  VALUES ?award {
    wd:Q255032 # Hugo Award
    wd:Q2576795 # Locus Award
    wd:Q266012 # Nebula Award
    wd:Q898527 # World Fantasy Award
    wd:Q1341487 # Philip K. Dick Award
    wd:Q1030402 # John W. Campbell Memorial Award
    wd:Q39060754 # Dragon Award
    wd:Q1666525 # International Fantasy Award
  }.
  ?item wdt:P50 ?author.
  ?item wdt:P136 ?genre.
  BIND(YEAR(?date) AS ?year)
  FILTER(?year <= 2017)

  ?item wdt:P7937 wd:Q8261. 
  
  # Check for identifiers for the work
  OPTIONAL { ?item wdt:P8383 ?work_id } # goodreads
  OPTIONAL { ?item wdt:P648 ?ol_id } # openlibrary
  OPTIONAL { ?item wdt:P1274 ?isfdb_id } # isfdb

  # Get text representations of the title, author, genre, and award 
  SERVICE wikibase:label {
    bd:serviceParam wikibase:language "[AUTO_LANGUAGE],mul,en".
    ?item rdfs:label ?title.
    ?author rdfs:label ?authorLabel.
    ?genre rdfs:label ?genreLabel.
    ?award rdfs:label ?awardLabel.
  }
}
GROUP BY ?item ?title ?year ?awardLabel ?status
ORDER BY DESC(?year)
"""


def get_results(endpoint_url, query):
    user_agent = "scifi-fantasy-awards (dev@zmraines.com) python/%s.%s" % (
        sys.version_info[0],
        sys.version_info[1],
    )

    headers = {"User-Agent": user_agent, "Accept": "application/sparql-results+json"}

    query = "\n".join(
        line.split("#")[0] for line in query.split("\n") if line.split("#")[0].strip()
    )
    params = {"query": query, "format": "json"}

    with httpx.Client() as client:
        response = client.get(
            endpoint_url, headers=headers, params=params, timeout=15.0
        )
        response.raise_for_status()
        return response.json()


results = get_results(endpoint_url, query)

with open("data/raw/wikidata_awards.json", "w") as f:
    for result in results["results"]["bindings"]:
        collapsed_result = {}
        for key, value in result.items():
            collapsed_result[key] = value["value"]

        print(json.dumps(collapsed_result), file=f)
