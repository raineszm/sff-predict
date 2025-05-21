-- Load the works with at least ${RATINGS_THRESHOLD} ratings
CREATE TABLE works AS
SELECT -- explicitly choose the columns we want to read
    -- and cast them to the correct types
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
    '${INPUT_WORKS}'
WHERE
    TRY_CAST(ratings_count as BIGINT) >= ${RATINGS_THRESHOLD}; -- only include works with at least ${RATINGS_THRESHOLD} ratings

-- Load the books that are part of the works we have kept
CREATE TABLE books AS
SELECT -- explicitly choose the columns we want to read
    -- and cast them to the correct types
    CAST(work_id as BIGINT) as work_id,
    CAST(book_id as BIGINT) as book_id,
    isbn,
    isbn13,
    asin,
    -- kindle_asin,
    title,
    title_without_series,
    -- TRY_CAST(text_reviews_count as BIGINT) as text_reviews_count,
    -- TRY_CAST(ratings_count as BIGINT) as ratings_count,
    -- TRY_CAST(average_rating as DOUBLE) as average_rating,
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
    -- url,
    -- link,
    -- image_url,
    series,
    CAST(popular_shelves as STRUCT(count BIGINT, name VARCHAR)[]) as popular_shelves,
    list_transform(authors, s -> s.author_id) as author_ids,
    -- CAST(similar_books as BIGINT[]) as similar_books,
FROM
    '${INPUT_BOOKS}'
SEMI JOIN works USING (work_id) -- only include books that are part of a work we have kept
WHERE
    language_code LIKE 'en%' -- only include English books
    AND work_id != '' -- Drop entries with no work ID since we don't know what work they belong to
;

-- Create a view of all tags for each work
-- and filter for only those works which are tagged as scifi/fantasy
-- enough times
CREATE VIEW tags_by_work_id AS
SELECT
    work_id,
    LIST (tags) as tags -- combine the tags into a list
FROM
    (
        -- this subsquery explodes each book into multiple rows,
        -- one for each popular shelf it is on
        SELECT
            work_id,
            UNNEST (popular_shelves) as tags
        FROM
            books
    )
GROUP BY
    work_id -- we then group by work_id to get the list of tags for each work
HAVING
    SUM(tags.count) FILTER ( -- filter for works with at least ${TAG_THRESHOLD} tags
        WHERE
            tags.name ~ '(?i)fantasy|sci(ence)?-fi'
    ) >= ${TAG_THRESHOLD};
-- and finally, see the SELECT clause above where we collapse the multiple rows,
-- back down to one row, with a list of tags

-- Add the tags to the works table
CREATE TABLE tagged_works AS
SELECT
    *
FROM
    works
INNER JOIN tags_by_work_id USING (work_id);

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

CREATE TABLE authors_by_work_id AS
WITH all_authors AS (
    SELECT 
        author_id,
        name as author_name
    FROM '${INPUT_AUTHORS}'
),
used_author_ids AS (
    SELECT
        work_id,
        UNNEST(author_ids) as author_id
    FROM combined_works
)
SELECT
    work_id, array_to_string(
        list_sort(LIST(author_name))
        , ', ') as authors
FROM
    used_author_ids
INNER JOIN all_authors USING (author_id)
GROUP BY work_id;

-- And save the augmented works to a parquet file
COPY (
    SELECT
        * EXCLUDE(author_ids)
    FROM
        combined_works
    INNER JOIN authors_by_work_id USING (work_id)
) TO '${OUTPUT_AUGMENTED_WORKS}' (FORMAT PARQUET, COMPRESSION ZSTD);