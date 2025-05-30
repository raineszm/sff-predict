import os
import pandas as pd
from typing import Optional
from googleapiclient.discovery import build
from tqdm.auto import tqdm


def _get_volume(client, query: Optional[str] = None) -> dict:
    return (
        client.volumes()
        .list(
            q=query,
            projection="LITE",
            printType="BOOKS",
            fields="totalItems,items(volumeInfo/description, volumeInfo/publishedDate,volumeInfo/industryIdentifiers)",
        )
        .execute()
    )


def lookup_volume(client, row: pd.Series) -> dict:
    if pd.notna(row.get("isbn13")) and row["isbn13"] != "":
        return _get_volume(client, f"isbn:{row['isbn13']}")
    elif pd.notna(row.get("isbn")) and row["isbn"] != "":
        return _get_volume(client, f"isbn:{row['isbn']}")
    else:
        first_author = row["authors"].split(",")[0]
        return _get_volume(client, f"intitle:{row['title']}+inauthor:{first_author}")


def get_volume_data(client, row: pd.Series) -> Optional[dict]:
    volume = lookup_volume(client, row)
    if volume["totalItems"] == 0:
        return None

    info = volume["items"][0]["volumeInfo"]

    if (
        "description" not in info
        or "publishedDate" not in info
        or "industryIdentifiers" not in info
    ):
        return None

    return {
        "description": info["description"],
        "publication_year": info["publishedDate"][:4],
        "isbn": info["industryIdentifiers"][0]["identifier"],
    }


def fill_incomplete_works(
    df: pd.DataFrame, api_key: Optional[str] = None
) -> pd.DataFrame:
    """Fill in missing data for incomplete works using the Google Books API."""

    # Identify rows with missing data
    incomplete_mask = (
        (df.description == "")
        | df.description.isna()
        | (df.publication_year.isna() & df.original_publication_year.isna())
    )

    # Fill missing data using Google Books API

    to_drop = []
    with build("books", "v1", developerKey=api_key) as client:
        for index in tqdm(
            df[incomplete_mask].index,
            total=incomplete_mask.sum(),
            desc="Filling missing data",
        ):
            volume_data = get_volume_data(client, df.loc[index])
            if volume_data is None:
                to_drop.append(index)
                continue

            df.loc[index, "description"] = volume_data["description"]
            df.loc[index, "publication_year"] = volume_data["publication_year"]
            df.loc[index, "isbn"] = volume_data["isbn"]

    df = df.drop(to_drop, axis=0)

    return df


if __name__ == "__main__":
    data = pd.read_parquet("data/sff_works_augmented.parquet")
    filled_data = fill_incomplete_works(data, os.environ.get("GOOGLE_API_KEY"))
    filled_data.to_parquet("data/sff_works_augmented_filled.parquet")
