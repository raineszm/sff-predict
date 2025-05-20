-- Load the works with at least ${RATINGS_THRESHOLD} ratings
CREATE TABLE works AS
SELECT
    CAST(work_id as BIGINT) as work_id,
    CAST(best_book_id as BIGINT) as best_book_id,
    -- TRY_CAST(books_count as BIGINT) as books_count,
    TRY_CAST(reviews_count as BIGINT) as reviews_count,
    TRY_CAST(text_reviews_count as BIGINT) as text_reviews_count,
    TRY_CAST(ratings_count as BIGINT) as ratings_count,
    -- TRY_CAST(ratings_sum as BIGINT) as ratings_sum,
    TRY_CAST(original_publication_year as INTEGER) as original_publication_year,
    TRY_CAST(original_publication_month as INTEGER) as original_publication_month,
    TRY_CAST(original_publication_day as INTEGER) as original_publication_day,
    original_title
    -- original_language_id,
    -- default_description_language_code,
    -- default_chaptering_book_id,
    -- media_type,
    -- rating_dist,
FROM
    read_ndjson_auto('${INPUT_WORKS}', convert_strings_to_integers = true)
WHERE
    ratings_count >= ${RATINGS_THRESHOLD};

ALTER TABLE works
ADD PRIMARY KEY (work_id);

CREATE TABLE books AS
SELECT
    CAST(work_id as BIGINT) as work_id,
    CAST(book_id as BIGINT) as book_id,
    isbn,
    isbn13,
    asin,
    -- kindle_asin,
    title,
    title_without_series,
    TRY_CAST(text_reviews_count as BIGINT) as text_reviews_count,
    TRY_CAST(ratings_count as BIGINT) as ratings_count,
    TRY_CAST(average_rating as DOUBLE) as average_rating,
    TRY_CAST(num_pages as BIGINT) as num_pages,
    -- format,
    publisher,
    TRY_CAST(publication_year as INTEGER) as publication_year,
    TRY_CAST(publication_month as INTEGER) as publication_month,
    TRY_CAST(publication_day as INTEGER) as publication_day,
    -- edition_information,
    language_code,
    country_code,
    CAST(is_ebook as BOOLEAN) as is_ebook,
    description,
    url,
    -- link,
    -- image_url,
    series,
    CAST(popular_shelves as STRUCT(count BIGINT, name VARCHAR)[]) as popular_shelves,
    CAST(authors as STRUCT(author_id BIGINT, role VARCHAR)[]) as authors,
    CAST(similar_books as BIGINT[]) as similar_books,
FROM
    '${INPUT_BOOKS}'
WHERE
    language_code LIKE 'en%'
    AND work_id != '' -- Drop entries with no work ID since we don't know what work they belong to
    AND CAST(work_id as BIGINT) IN (
        SELECT
            work_id
        FROM
            works
    );


ALTER TABLE books
ADD PRIMARY KEY (book_id);

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
    AND SUM(tags.count) FILTER (
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

ALTER TABLE tagged_works
ADD PRIMARY KEY (work_id);

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
    * EXCEPT (books.work_id)
FROM
    tagged_works
INNER JOIN books ON tagged_works.best_book_id = books.book_id;

ALTER TABLE combined_works
ADD PRIMARY KEY (work_id);

-- And save the augmented works to a parquet file
COPY (
    SELECT
        *
    FROM
        combined_works
) TO '${OUTPUT_AUGMENTED_WORKS}' (FORMAT PARQUET, COMPRESSION ZSTD);