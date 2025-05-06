#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Sample top N rows from a large edgelist and interest CSV."
    )
    parser.add_argument("n", type=int, help="Number of data rows to extract")
    parser.add_argument(
        "--edgelist", "-e", default="full.edgelist",
        help="Path to the full edgelist file (space-delimited, no header)"
    )
    parser.add_argument(
        "--interest", "-i", default="full_interests.csv",
        help="Path to the full interests CSV (comma-delimited, with header)"
    )
    parser.add_argument(
        "--out-edge", "-o", default="sample.edgelist",
        help="Path for the sampled edgelist output"
    )
    parser.add_argument(
        "--out-interest", "-r", default="sample_interests.csv",
        help="Path for the sampled interests output"
    )
    args = parser.parse_args()

    # Sample edgelist (no header)
    with open(args.edgelist, 'r') as fin, open(args.out_edge, 'w') as fout:
        for i, line in enumerate(fin):
            if i >= args.n:
                break
            fout.write(line)

    # Sample interests CSV (assume first line is header)
    with open(args.interest, 'r') as fin, open(args.out_interest, 'w') as fout:
        header = next(fin)
        fout.write(header)
        for i, line in enumerate(fin):
            if i >= args.n:
                break
            fout.write(line)

if __name__ == "__main__":
    main()
