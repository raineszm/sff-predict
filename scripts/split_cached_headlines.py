import pandas as pd


df = pd.read_json("data/raw/headlines.json", lines=True)

for (year, month), headlines in df.groupby(["year", "month"]):
    headlines.to_json(
        f"data/raw/headlines/{year}-{month}.json", orient="records", lines=True
    )
