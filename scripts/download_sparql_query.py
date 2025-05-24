from snakemake.script import snakemake
import sys
import httpx
import json

endpoint_url = "https://query.wikidata.org/sparql"


def get_results(endpoint_url, query):
    user_agent = "scifi-fantasy-awards (dev@zmraines.com) python/%s.%s" % (
        sys.version_info[0],
        sys.version_info[1],
    )

    headers = {"User-Agent": user_agent, "Accept": "application/sparql-results+json"}

    params = {"query": query, "format": "json"}

    with httpx.Client() as client:
        response = client.get(
            endpoint_url, headers=headers, params=params, timeout=15.0
        )
        response.raise_for_status()
        return response.json()


with open("scripts/queries/wikidata_awards.sparql", "r") as f:
    query = "\n".join(line.split("#")[0] for line in f if line.split("#")[0].strip())


def download_query(query_path, output_path):
    with open(query_path, "r") as f:
        query = "\n".join(
            line.split("#")[0] for line in f if line.split("#")[0].strip()
        )

    results = get_results(endpoint_url, query)

    with open(output_path, "w") as f:
        for result in results["results"]["bindings"]:
            collapsed_result = {}
            for key, value in result.items():
                collapsed_result[key] = value["value"]
            print(json.dumps(collapsed_result), file=f)


download_query(snakemake.input.sparql, snakemake.output.json)
