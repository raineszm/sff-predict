-- Load the works with at least ${RATINGS_THRESHOLD} ratings
CREATE TABLE works AS
SELECT
    *
FROM
    '${INPUT_WORKS}'
WHERE
    CAST(ratings_count AS INTEGER) >= ${RATINGS_THRESHOLD};

CREATE TABLE books AS
SELECT
    *
FROM
    '${INPUT_BOOKS}'
WHERE
    language_code LIKE 'en%'
    AND work_id IN (
        SELECT
            work_id
        FROM
            works
    );



-- Create a table of all tags for each work
-- and filter for only those which are tagged as scifi/fantasy
-- enough times
CREATE TABLE work_tags AS
SELECT
    work_id,
    LIST (tags)
FROM
    (
        SELECT
            work_id,
            UNNEST (popular_shelves) as tags
        FROM
            books
    )
GROUP BY
    work_id
HAVING
    work_id IN (
        SELECT
            work_id
        FROM
            works
    )
    AND SUM(CAST(tags.count as INTEGER)) FILTER (
        WHERE
            tags.name ~ '(?i)fantasy|sci(ence)?-fi'
    ) >= ${TAG_THRESHOLD};

-- Combine the selected works with 
-- the date in the works table
CREATE TABLE tagged_works AS
SELECT
    *
FROM
    work_tags
INNER JOIN works USING (work_id);

-- Save the filtered works to a parquet file
COPY (
    SELECT
        *
    FROM
        tagged_works
) TO '${OUTPUT_SELECTED_WORKS}' (FORMAT PARQUET, COMPRESSION ZSTD);

-- Augment the works data with data from a representative edition
-- of the work
CREATE TABLE combined_works AS
SELECT
    *
FROM
    tagged_works
INNER JOIN books ON tagged_works.best_book_id = books.book_id;

-- And save the augmented works to a parquet file
COPY (
    SELECT
        *
    FROM
        combined_works
) TO '${OUTPUT_AUGMENTED_WORKS}' (FORMAT PARQUET, COMPRESSION ZSTD);