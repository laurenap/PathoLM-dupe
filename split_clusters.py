import csv, random
from pathlib import Path
random.seed(1337)

# Load header -> label (built earlier)
lab = {}
with open("manifests/header_label.tsv") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        lab[row["seq_id"]] = int(row["label"])

def split_from_clusters(tsv_path: str, outpref: str):
    # mmseqs createtsv format: <rep_seq_id>\t<member_seq_id>
    clusters = {}
    with open(tsv_path) as f:
        for rep, mem in csv.reader(f, delimiter="\t"):
            clusters.setdefault(rep, []).append(mem)

    # Majority label per cluster (only members we have labels for)
    items = []
    for cid, members in clusters.items():
        labels = [lab[m] for m in members if m in lab]
        if not labels:
            continue
        maj = 1 if sum(labels) >= len(labels)/2 else 0
        items.append((cid, maj, members))

    pos = [x for x in items if x[1] == 1]
    neg = [x for x in items if x[1] == 0]
    random.shuffle(pos); random.shuffle(neg)

    def cut(gs, frac):
        n = int(len(gs) * frac)
        return gs[:n], gs[n:]

    pos_tr, pos_te = cut(pos, 0.8)
    neg_tr, neg_te = cut(neg, 0.8)

    train = pos_tr + neg_tr
    test  = pos_te + neg_te
    random.shuffle(train); random.shuffle(test)

    # 10% of train -> val
    nval = max(1, int(0.1 * len(train)))
    val, train = train[:nval], train[nval:]

    def dump(name, rows):
        out = Path(f"splits/{outpref}_{name}.txt")
        with out.open("w") as f:
            for _, _, members in rows:
                for m in members:
                    f.write(m + "\n")

    dump("train", train); dump("val", val); dump("test", test)
    print(outpref, "clusters:", len(items),
          "-> train:", len(train), "val:", len(val), "test:", len(test))

for t in (80, 60, 40):
    split_from_clusters(f"splits/clusters_{t}.tsv", f"split{t}")
