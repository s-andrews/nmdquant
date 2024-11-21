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

    write_output(splices, options.bam, options.outfile, options.nounmeasured)

    print("Run complete")

def write_output(splices, files, outfile, ignoreunmeasured):
    print("Writing output to",outfile)
    with open(outfile,"wt",encoding="utf8") as out:
        headers = ["INTRON","DIRECTION","GENE","NMD"]
        headers.extend(files)

        print("\t".join(headers), file=out)

        for intron in splices.keys():

            # If they've set --nounmeasured then we don't print lines
            # where the quantitations are all zero
            if ignoreunmeasured:
                reject_line = True
                for value in splices[intron]["quantitations"]:
                    if value > 0:
                        reject_line = False
                        break
                if reject_line:
                    continue

            line = [intron,splices[intron]["direction"],splices[intron]["gene"],splices[intron]["nmd"]]
            line.extend(splices[intron]["quantitations"])

            print("\t".join([str(x) for x in line]), file=out)


def quantitate_bam(bam, splices):
    print("Quantitating",bam)
    # We first need to add a zero value to every splice in the dataset
    # then we can parse through the file and increment every time we
    # find a splice which is listed in the file.
    for splice in splices:
        splices[splice]["quantitations"].append(0)

    inbam = pysam.AlignmentFile(bam,"rb")

    for read_count,read in enumerate(inbam.fetch(until_eof=True)):
        if read_count % 1000000 == 0:
            print("Processed",int(read_count/1000000),"million reads")

        # We don't want non unique mappings
        if read.mapping_quality < 20:
            continue

        splice_lengths = []
        for operation,operation_length in read.cigartuples:
            if operation == 3:
                splice_lengths.append(operation_length)

        if not splice_lengths:
            continue

        blocks = read.blocks
        for i in range(1,len(blocks)):
            dist = blocks[i][0] - blocks[i-1][1]
            if dist == splice_lengths[0]:
                splice_lengths.pop(0)
                intron = f"{read.reference_name}:{blocks[i-1][1]+1}-{blocks[i][0]}"
                if intron in splices:
                    splices[intron]["quantitations"][-1] += 1

                if not splice_lengths:
                    break
                

def read_splice_sites(file):
    print("Extracting features from",file)
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

    last_chr = ""

    for line in infh:
        line = line.strip()

        if line.startswith("#"):
            continue

        sections = line.split("\t")

        if sections[0] != last_chr:
            print("Reading features from chr",sections[0])
            last_chr = sections[0]

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
            gene_name=annotations["gene_id"]
            if "gene_name" in annotations:
                gene_name = annotations["gene_name"]
            transcripts[annotations["transcript_id"]] = {
                "chromosome": sections[0],
                "direction" : sections[6],
                "gene" : gene_name,
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

    for transcript_id in transcripts:
        transcript_exons = transcripts[transcript_id]["exons"]
        cds_exons = cds[transcript_id]

        if not cds_exons:
            continue

        cds_start = None
        cds_end = None

        # We need to put the exons in order.  This is based on the
        # orientation of the gene

        if transcripts[transcript_id]["direction"] == "+":
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
            if transcripts[transcript_id]['direction'] == "+":
                intron = f"{transcripts[transcript_id]['chromosome']}:{transcript_exons[i-1][1]+1}-{transcript_exons[i][0]-1}"
                if cds_end is None or transcript_exons[i-1][1]+1 > cds_end:
                    is_nmd = True

            else:
                intron = f"{transcripts[transcript_id]['chromosome']}:{transcript_exons[i][1]+1}-{transcript_exons[i-1][0]-1}"
                if cds_end is None or transcript_exons[i][1]-1 < cds_end:
                    is_nmd = True

            if not intron in splice_sites:
                splice_sites[intron] = {"nmd":is_nmd, "direction":transcripts[transcript_id]['direction'], "quantitations":[], "gene": transcripts[transcript_id]['gene']}

            if not is_nmd:
                splice_sites[intron]["nmd"] = False        

    return splice_sites

def get_options():
    parser = argparse.ArgumentParser(description="A program to quantify nonsense mediated decay in RNA-Seq data")

    parser.add_argument("gtf", help="GTF file of annotations")
    parser.add_argument("bam", nargs="+", help="BAM files to quantitate")
    parser.add_argument("--outfile",required=True, help="Output file name")
    parser.add_argument("--nounmeasured", default=False, action="store_true", help="Don't report introns with no counts in any sample")

    return parser.parse_args()


if __name__ == "__main__":
    main()