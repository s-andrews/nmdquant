Quantitation of Nonsense Mediated Decay
=======================================

[Nonsense Mediated Decay](https://en.wikipedia.org/wiki/Nonsense-mediated_decay), or NMD for short, is a surveilance mechanims in cells which identifies and removes protein coding RNA molecules which contain a premature termination codon.  Any transcript which contains splice junctions, and where the coding sequence within the transcript does not extend beyond the final splice junction is susceptible to degradation via this mechanism.

![NMD illustration](https://upload.wikimedia.org/wikipedia/commons/5/52/NMD_-_Nonsense-mediated_decay.png)

Under some circumstances this process may be inhibited, leading to an increase in the number of non-productive transcripts within the mRNA population.

The ```nmdquant``` program aims to quantitate the level of nonsense mediated decay target transcripts within a set of RNA-Seq samples.  It takes in BAM files of RNA-Seq data which have been mapped using a splicing aware mapper such as [Hisat2](http://daehwankimlab.github.io/hisat2/) or [STAR](https://github.com/alexdobin/STAR) and produces counts of reads crossing all splice boundaries in a supplied GTF annotation file.  The splice boundaries are annotated as to whether they are exclusively found in NMD target transcripts.

From the output of ```nmdquant``` you can visualise and analyse the levels of NMD in different samples.

Installation
------------

```nmdquant``` is a python program.  To install this on your system you can download the code from this repository.  You will need to use any recent version of python (3.7 onwards), and you will need to install the ```pysam``` package using ```python3 -m pip install pysam```.  After that you can run the ```nmdqant.py``` script to start the analysis.


Usage
-----

```
$ ./nmdquant.py --help
usage: nmdquant.py [-h] --outfile OUTFILE [--nounmeasured] gtf bam [bam ...]

A program to quantify nonsense mediated decay in RNA-Seq data

positional arguments:
  gtf                GTF file of annotations
  bam                BAM files to quantitate

options:
  -h, --help         show this help message and exit
  --outfile OUTFILE  Output file name
  --nounmeasured     Don't report introns with no counts in any sample
```

Inputs to the program are:

1. One or more BAM files of RNA-Seq data mapped with a splicing aware read mapper to a suitable reference genome.  The BAM files do not need to be sorted.
2. A GTF annotation file for the genome assembly against which the reads were aligned.  The program has been tested with GTF files coming from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) but other sources should also work.  The GTF file can optionall be gzip compressed

The output of the program is a tab delimited text file of the raw counts for each BAM file of the number of reads showing each of the annotated splice sites in the GTF file.  Columns in the output are:

1. INTRON (location of the intron, eg 1:110012904-110026000)  
2. DIRECTION (the strand from which the intron came + or -). Quantitation is not directional
3. GENE (the name of the gene from which the intron comes, or the gene id if no name is provided)
4. NMD (True or False, whether this intron is solely found in nmd susceptible transcripts)
5. BAM Names - the remaining columns represent the names of the BAM files being quantitated


Method
------

The program initially goes through the GTF file assembling the exons for each transcript and CDS feature. It then extracts the introns between the annotated exons, and tags these as NMD if they occur after the end of the CDS in a given transcript.  If an intron is ambiguous (NMD target in one transcript, but not in another), then it is **not** flagged as NMD.

After assmbling the list of annotated introns from the GTF the program then iterates through the BAM files, identifying all observed splice junctions and relating them to the annotated introns from the GTF file.  Matches between the GTF and the BAM files are then counted.  Any splices observed in the BAM files but not annotated in the GTF are ignored.

Finally the results are written to a tab delimited text file.

Downstream Analysis
-------------------

The end point of the program is the table of counts for each annotated intron.  Downstream analysis of these values can be performed in any suitable environment.  One example of such an analysis in R is given in the ```analysis``` subfolder of this repository.

[View Example Analysis](https://htmlpreview.github.io/?https://raw.githubusercontent.com/s-andrews/nmdquant/refs/heads/main/analysis/nmd_quantitation_analysis.html)

