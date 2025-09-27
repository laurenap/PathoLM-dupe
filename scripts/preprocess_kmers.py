# preprocess_kmers.py
from pathlib import Path
from textwrap import shorten

# =========================
# CONFIG
# =========================
K = 6  # k-mer length
BASE = Path(__file__).resolve().parent
DATA = BASE / "data_collection"  # <<— uses data_collection (no space)

# Expected directory layout (under data_collection/)
DIRS = {
    "bact_pathogenic": DATA / "Bacteria_ESKAPEE" / "pathogenic set",
    "bact_nonpath":    DATA / "Bacteria_ESKAPEE" / "non-pathogenic set",
    "virus_path":      DATA / "Viruses" / "Human-pathogenic",
    "virus_nonpath":   DATA / "Viruses" / "Non-pathogenic comparators",
}

OUT = BASE / "tokenized"
OUT.mkdir(exist_ok=True)
PLAIN_CORPUS = OUT / "corpus.txt"
KMER_CORPUS  = OUT / f"corpus_k{K}.txt"

# =========================
# HELPERS
# =========================
def fasta_iter(p: Path):
    """Yield (header, seq) tuples from a fasta/fna/fa file (plaintext)."""
    header, seq = None, []
    with p.open() as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)

def clean_seq(s: str) -> str:
    s = s.upper()
    return "".join(c if c in "ACGTN" else "N" for c in s)

def kmers(s: str, k: int):
    n = len(s)
    for i in range(0, max(0, n - k + 1)):
        yield s[i:i+k]

def find_genomes(root: Path):
    """Recursively find fasta-like files under root."""
    # Only plaintext files; add gzip handling later if needed
    valid = {".fa", ".fna", ".fasta"}
    return [p for p in root.rglob("*") if p.suffix.lower() in valid]

def exists_or_try_variants(path: Path) -> Path | None:
    """
    Return a path that exists, trying common filename gotchas:
    - trailing spaces
    - non-breaking spaces
    - en dash vs hyphen
    """
    if path.exists():
        return path
    # Try NBSP and trailing-space normalization on each component
    parts = list(path.parts)
    base = Path(parts[0])
    for comp in parts[1:]:
        # search within base for a child that loosely matches comp
        candidates = []
        if base.exists():
            for child in base.iterdir():
                if not child.exists():
                    continue
                norm_child = child.name.replace("\u00A0", " ").replace("–", "-").rstrip()
                norm_comp  = comp.replace("\u00A0", " ").replace("–", "-").rstrip()
                if norm_child.lower() == norm_comp.lower():
                    candidates.append(child)
        base = candidates[0] if candidates else base / comp
    return base if base.exists() else None

# =========================
# MAIN
# =========================
def main():
    print(f"[info] Working dir: {BASE}")
    print(f"[info] Data dir:     {DATA}")

    totals = {}
    files_map = {}

    for key, raw_path in DIRS.items():
        # resolve path, trying to heal common naming issues
        path = raw_path if raw_path.exists() else exists_or_try_variants(raw_path)
        if not path or not path.exists():
            print(f"[warn] Missing folder: {raw_path}")
            totals[key] = 0
            files_map[key] = []
            continue

        files = find_genomes(path)
        files_map[key] = files
        totals[key] = len(files)

        status = "[ok]" if files else "[warn]"
        print(f"{status} {key}: {path}  -> {len(files)} fasta files")
        for p in files[:5]:
            print(f"      - {p}")

    if all(v == 0 for v in totals.values()):
        print("[fatal] Found zero fasta files. Check your folder names and location of fasta files.")
        return

    # Build corpora
    plain_lines = []
    kmer_lines  = []

    for label, files in files_map.items():
        # Binary label: 1 = pathogenic, 0 = non-pathogenic
        y = 1 if ("path" in label and "non" not in label) else 0
        for fp in files:
            try:
                for hdr, seq in fasta_iter(fp):
                    s = clean_seq(seq)
                    if not s:
                        continue
                    # plain sequence corpus
                    plain_lines.append(f"{y}\t{shorten(hdr, width=80, placeholder='…')}\t{s}\n")
                    # k-mer corpus
                    kline = " ".join(kmers(s, K))
                    if kline:
                        kmer_lines.append(f"{y}\t{shorten(hdr, width=80, placeholder='…')}\t{kline}\n")
            except Exception as e:
                print(f"[err] Failed on {fp}: {e}")

    PLAIN_CORPUS.write_text("".join(plain_lines))
    KMER_CORPUS.write_text("".join(kmer_lines))

    print("\n[done] Preprocessing complete!")
    print(f"       Wrote plain corpus -> {PLAIN_CORPUS}  ({len(plain_lines)} sequences)")
    print(f"       Wrote k-mer corpus -> {KMER_CORPUS}   ({len(kmer_lines)} sequences)")
    print("       Tip: open the files and verify the first few lines aren’t empty.")

if __name__ == "__main__":
    main()