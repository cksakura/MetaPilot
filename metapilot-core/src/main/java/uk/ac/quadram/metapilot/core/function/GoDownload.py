from pathlib import Path

obo = Path("go-basic.obo")
out = Path("GO.tsv")

records = []
current_id = None
current_name = None

with obo.open(encoding="utf-8") as f:
    for line in f:
        line = line.rstrip("\n")

        if line == "[Term]":
            current_id = None
            current_name = None

        elif line.startswith("id: GO:"):
            current_id = line.replace("id: ", "")

        elif line.startswith("name: "):
            current_name = line.replace("name: ", "")

        elif line == "":
            if current_id and current_name:
                records.append((current_id, current_name))
            current_id = None
            current_name = None

with out.open("w", encoding="utf-8") as w:
    for go_id, name in records:
        w.write(f"{go_id}\t{name}\n")