from snakemake.script import snakemake
import duckdb

duckdb.sql(
    f"CREATE OR REPLACE TABLE openlibrary_ids AS FROM '{snakemake.input['openlibrary_ids']}'"
)

descriptions = duckdb.sql(
    """
SELECT DISTINCT openlibrary_ids.work_qid, openlibrary_works.description
FROM openlibrary_ids
LEFT JOIN '{ol_works}' as openlibrary_works
ON '/works/' || openlibrary_ids.openlibrary_id = openlibrary_works.key
WHERE openlibrary_works.description IS NOT NULL
""".format(
        ol_works=snakemake.input["ol_works"],
    )
)


n_works = duckdb.sql(
    "SELECT DISTINCT COUNT(work_qid) FROM '{nominated_novels}';".format(
        nominated_novels=snakemake.input["nominated_novels"]
    )
).fetchone()[0]
print("Descriptions found for {}/{} works".format(len(descriptions), n_works))

descriptions.write_csv(snakemake.output["descriptions"])
