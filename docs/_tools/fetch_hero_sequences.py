"""Refresh the protein sequences that drift behind the landing page hero.

The docs build never touches the network, so the sequences are stored in
docs/_static/landing/sequences.json. Run this script by hand to regenerate
that file, for example after changing ENTRIES:

    python docs/_tools/fetch_hero_sequences.py

Each record keeps the PDB entry, the chain label and the source URL, so
anyone can trace a sequence on the page back to its structure.
"""

import json
import urllib.request
from pathlib import Path

# Public structures from the lab's research (see docs/design/PRODUCT.md).
ENTRIES = ["7MU9", "4QTQ", "7LKQ", "4FOU", "8PX4"]

FASTA_URL = "https://www.rcsb.org/fasta/entry/{entry}"
OUT = Path(__file__).resolve().parents[1] / "_static" / "landing" / "sequences.json"


def fetch_fasta(entry):
    with urllib.request.urlopen(FASTA_URL.format(entry=entry), timeout=30) as resp:
        return resp.read().decode("utf-8")


def parse_fasta(text):
    """Yield (header, sequence) pairs."""
    header, seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def records():
    seen = set()
    for entry in ENTRIES:
        for header, seq in parse_fasta(fetch_fasta(entry)):
            if seq in seen:  # the same chain appears in more than one entry
                continue
            seen.add(seq)
            fields = header.split("|")
            yield {
                "entry": entry,
                "chain": fields[1] if len(fields) > 1 else "",
                "molecule": fields[2] if len(fields) > 2 else "",
                "organism": fields[3] if len(fields) > 3 else "",
                "source": FASTA_URL.format(entry=entry),
                "sequence": seq,
            }


def main():
    data = list(records())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, indent=2) + "\n")
    print(f"wrote {len(data)} sequences to {OUT}")


if __name__ == "__main__":
    main()
