import requests
import pandas as pd

url = "http://bigg.ucsd.edu/api/v2/universal/reactions"
data = requests.get(url).json()

rows = []
for r in data["results"]:
    rows.append({
        "BiGG_Reaction": r["bigg_id"],
        "Description": r["name"]
    })

df = pd.DataFrame(rows)
df.to_csv("BiGG_Reaction.tsv", sep="\t", index=False)