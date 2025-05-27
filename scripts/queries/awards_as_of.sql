CREATE OR REPLACE TABLE cumulative_awards AS WITH awards_per_year AS (
        SELECT author_qid,
            authorLabel,
            year,
            COUNT(awardLabel) as award_count
        FROM '${INPUT_AWARDS}'
        GROUP BY author_qid,
            authorLabel,
            year
    )
SELECT author_qid,
    authorLabel,
    year,
    SUM(award_count) OVER (
        PARTITION BY author_qid
        ORDER BY year
    ) as awards_as_of_year
FROM awards_per_year
ORDER BY authorLabel,
    year;
COPY (
    FROM cumulative_awards
) TO '${OUTPUT_CUMULATIVE}';