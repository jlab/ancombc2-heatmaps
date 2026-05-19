from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Sequence


RANK_TO_PREFIX = {
    "d": "d__",
    "k": "k__",
    "p": "p__",
    "c": "c__",
    "o": "o__",
    "f": "f__",
    "g": "g__",
}


def clean_tax_piece(x) -> Optional[str]:
    if x is None:
        return None
    x = str(x).strip()
    return x if x else None


def is_empty_tax(x) -> bool:
    if x is None:
        return True

    stripped = re.sub(r"^[a-z]__+", "", str(x)).strip().lower()
    return stripped in {
        "",
        "_",
        "uncultured",
        "unclassified",
        "unknown",
        "ambiguous_taxa",
    }


def normalize_taxon_label(raw_tax) -> str:
    """
    Convert QIIME-style taxonomy strings to:
    top_rank; family; genus

    Example:
    d__Bacteria; p__Verrucomicrobiota; ...; f__Akkermansiaceae; g__Akkermansia
    -> p__Verrucomicrobiota; f__Akkermansiaceae; g__Akkermansia
    """
    if raw_tax is None:
        return "_; _; _"

    if isinstance(raw_tax, (list, tuple)):
        parts = [clean_tax_piece(x) for x in raw_tax]
    else:
        parts = [clean_tax_piece(x) for x in str(raw_tax).split(";")]

    parts = [x for x in parts if x is not None]

    tax_map = {}
    for part in parts:
        part = part.strip()
        if "__" not in part:
            continue

        prefix = part.split("__", 1)[0].lower()
        if not is_empty_tax(part):
            tax_map[prefix] = part

    top_rank = "_"
    for rank in ["o", "c", "p", "k", "d"]:
        if rank in tax_map:
            top_rank = tax_map[rank]
            break

    family = tax_map.get("f", "_")
    genus = tax_map.get("g", "_")

    return f"{top_rank}; {family}; {genus}"


def parse_label(label: str) -> Dict[str, str]:
    parts = [x.strip() for x in str(label).split(";")]
    while len(parts) < 3:
        parts.append("_")

    return {
        "top": parts[0],
        "family": parts[1],
        "genus": parts[2],
    }


def normalize_query(query: str) -> Dict[str, str]:
    query = str(query).strip()

    if ";" in query:
        return {"mode": "exact", "value": query}

    match = re.match(r"([a-zA-Z])_(.+)", query)
    if not match:
        raise ValueError(
            "Invalid taxon query. Use e.g. 'g_Akkermansia', "
            "'f_Akkermansiaceae', 'p_Verrucomicrobiota', or a full exact label."
        )

    rank = match.group(1).lower()
    name = match.group(2)

    if rank not in RANK_TO_PREFIX:
        raise ValueError(f"Unsupported taxonomic rank prefix: {rank}")

    return {
        "mode": "rank",
        "rank": rank,
        "value": RANK_TO_PREFIX[rank] + name,
    }


def match_taxa(labels: Iterable[str], query: str) -> List[str]:
    spec = normalize_query(query)
    matches = []

    for label in labels:
        parsed = parse_label(label)

        if spec["mode"] == "exact":
            if label == spec["value"]:
                matches.append(label)
            continue

        rank = spec["rank"]
        value = spec["value"]

        if rank == "f" and parsed["family"] == value:
            matches.append(label)
        elif rank == "g" and parsed["genus"] == value:
            matches.append(label)
        elif rank in ["o", "c", "p", "k", "d"] and parsed["top"] == value:
            matches.append(label)

    return matches


def list_available_queries(labels: Sequence[str]) -> Dict[str, List[str]]:
    families = set()
    genera = set()
    phyla = set()

    for label in labels:
        parsed = parse_label(label)

        family = parsed["family"]
        genus = parsed["genus"]
        top = parsed["top"]

        if family.startswith("f__") and not is_empty_tax(family):
            families.add("f_" + re.sub(r"^f__", "", family))

        if genus.startswith("g__") and not is_empty_tax(genus):
            genera.add("g_" + re.sub(r"^g__", "", genus))

        if top.startswith("p__") and not is_empty_tax(top):
            phyla.add("p_" + re.sub(r"^p__", "", top))

    return {
        "family_queries": sorted(families),
        "genus_queries": sorted(genera),
        "phylum_queries": sorted(phyla),
    }
