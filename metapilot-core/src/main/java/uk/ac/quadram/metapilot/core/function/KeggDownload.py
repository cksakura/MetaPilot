import pandas as pd

url = "https://rest.kegg.jp/list/ko"

ko = pd.read_csv(url, sep="\t", header=None, names=["KO", "Description"])
ko["KO"] = ko["KO"].str.replace("ko:", "", regex=False)

ko.to_csv("KEGG_KO_description.tsv", sep="\t", index=False)