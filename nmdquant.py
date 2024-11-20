#!/usr/bin/env python3
import argparse
import pysam
import gzip

options = None

def main():
    global options
    options = get_options()

    breakpoint()
    splices = read_splice_sites(options.gtf)


def read_splice_sites(file):
    infh = None
    splice_sites = {}
    if str(file).tolower().endswith(".gz"):
        infh = gzip.open(file, "rt", encoding="utf8")
    else:
        infh = open(file,"rt",encoding="utf8")



    infh.close()

    return splice_sites

def get_options():
    parser = argparse.ArgumentParser(description="A program to quantify nonsense mediated decay in RNA-Seq data")

    parser.add_argument("gtf", help="GTF file of annotations")
    parser.add_argument("bam", nargs="+", help="BAM files to quantitate")

    return parser.parse_args()


if __name__ == "__main__":
    main()