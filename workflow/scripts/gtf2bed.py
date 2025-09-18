#!/usr/bin/env python3
"""Convert a GTF annotation to BED12 format using gffutils."""

import argparse

import gffutils


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input GTF annotation file.")
    parser.add_argument("--bed", required=True, help="Output BED file path.")
    parser.add_argument("--db", required=True, help="Path to write the intermediate gffutils database.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    db = gffutils.create_db(
        args.input,
        dbfn=args.db,
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )

    with open(args.bed, "w") as outfileobj:
        for tx in db.features_of_type("transcript", order_by="start"):
            bed = [s.strip() for s in db.bed12(tx).split("\t")]
            bed[3] = tx.id
            outfileobj.write("{}\n".format("\t".join(bed)))


if __name__ == "__main__":
    main()
