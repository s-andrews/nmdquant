#!/usr/bin/env python3
import argparse
import pysam
import gzip

options = None

def main():
    global options
    options = get_options()

    splices = read_splice_sites(options.gtf)

    for bam in options.bam:
        quantitate_bam(bam,splices)

    write_output(splices, options.bam, options.outfile)

def write_output(splices, files, outfile):
    pass


def quantitate_bam(bam, splices):
    pass

def read_splice_sites(file):
    infh = None
    if str(file).lower().endswith(".gz"):
        infh = gzip.open(file, "rt", encoding="utf8")
    else:
        infh = open(file,"rt",encoding="utf8")

    # We want to identify all splice sites in the data as that's 
    # what we're going to quantitate.  We also want to flag those
    # which are only in NMD transcripts, that is transcripts where
    # the splice junction is after the end of any CDS which is in 
    # there.
        
    transcripts = {}
    cds = {}

    last_chr = "1"

    for line in infh:
        line = line.strip()

        if line.startswith("#"):
            continue

        sections = line.split("\t")

        if sections[0] != last_chr:
            print("Reading features from ",sections[0])
            last_chr = sections[0]
            ## TESTING ONLY
            break

        if not (sections[2] == "exon" or sections[2]=="CDS"):
            continue

        temp = sections[8].strip().split(";")
        annotations = {}
        for i in range(len(temp)):
            if not temp[i]:
                continue
            temp[i] = temp[i].strip().split(" ",1)
            annotations[temp[i][0]] = temp[i][1].replace('"','')


        if not annotations["transcript_id"] in transcripts:
            transcripts[annotations["transcript_id"]] = {
                "chromosome":sections[0],
                "direction" : sections[6],
                "exons" : []
            }
            cds[annotations["transcript_id"]] = []

        if sections[2] == "exon":
            transcripts[annotations["transcript_id"]]["exons"].append((int(sections[3]),int(sections[4])))

        elif sections[2] == "CDS":
            cds[annotations["transcript_id"]].append((int(sections[3]),int(sections[4])))

        else:
            raise Exception("Unknown feature type"+sections[2])

    infh.close()

    # Now we need to assemble the list of introns and flag them as nmd or valid
    # nmd exons are after the cds.  An intron which is nmd in one isoform might
    # not be in another so any ambiguous ones are treated as valid
    splice_sites = {}

    for transcipt_id in transcripts:
        transcript_exons = transcripts[transcipt_id]["exons"]
        cds_exons = cds[transcipt_id]

        if not cds_exons:
            continue

        cds_start = None
        cds_end = None

        # We need to put the exons in order.  This is based on the
        # orientation of the gene

        if transcripts[transcipt_id]["direction"] == "+":
            transcript_exons.sort(key=lambda x:x[0])
            for c in cds_exons:
                if cds_start  is None or c[0] < cds_start:
                    cds_start = c[0]
                if cds_end  is None or c[1] > cds_end:
                    cds_end = c[1]


        else:
            transcript_exons.sort(key=lambda x:x[0], reverse=True)
            for c in cds_exons:
                if cds_start  is None or c[1] > cds_start:
                    cds_start = c[1]
                if cds_end  is None or c[0] < cds_end:
                    cds_end = c[0]

        for i in range(1,len(transcript_exons)):
            is_nmd = False
            if transcripts[transcipt_id]['direction'] == "+":
                intron = f"{transcripts[transcipt_id]['chromosome']}:{transcript_exons[i-1][1]+1}-{transcript_exons[i][0]-1}:{transcripts[transcipt_id]['direction']}"
                if cds_end is None or transcript_exons[i-1][1]+1 > cds_end:
                    is_nmd = True

            else:
                intron = f"{transcripts[transcipt_id]['chromosome']}:{transcript_exons[i][1]+1}-{transcript_exons[i-1][0]-1}:{transcripts[transcipt_id]['direction']}"
                if cds_end is None or transcript_exons[i][1]-1 < cds_end:
                    is_nmd = True

            if not intron in splice_sites:
                splice_sites[intron] = {"nmd":is_nmd}

            if not is_nmd:
                splice_sites[intron]["nmd"] = False        

    return splice_sites

def get_options():
    parser = argparse.ArgumentParser(description="A program to quantify nonsense mediated decay in RNA-Seq data")

    parser.add_argument("gtf", help="GTF file of annotations")
    parser.add_argument("bam", nargs="+", help="BAM files to quantitate")
    parser.add_argument("--outfile",required=True, help="Output file name")

    return parser.parse_args()


if __name__ == "__main__":
    main()