#What sequencing depth is your alignment?

To calculates sequencing depth of coverage: 

C=(read length)(number of reads)/Genome length

Genome length (number of characters in the fasta file):

    ljcohen@farm:~/ECL243/Assignment2$ grep -v "^>" Zea_mays.AGPv3.30.dna_sm.chromosome.10.fa | wc -m
    152126075


Number of reads:

    ljcohen@farm:~/ECL243/Assignment2$ zcat Zmays_R1.fq.gz | wc -l
    1132416
    ljcohen@farm:~/ECL243/Assignment2$ zcat Zmex_R1.fq.gz | wc -l
    5346216
    ljcohen@farm:~/ECL243/Assignment2$ zcat Tdip_R1.fq.gz | wc -l
    17074344


Read length:

    ljcohen@farm:~/ECL243/Assignment2$ Zmex_length=$(zcat Zmex_R1.fq.gz | sed -n '2p')
    ljcohen@farm:~/ECL243/Assignment2$ echo ${#Zmex_length}
    76
    ljcohen@farm:~/ECL243/Assignment2$ Zmays_length=$(zcat Zmays_R1.fq.gz | sed -n '2p')
    ljcohen@farm:~/ECL243/Assignment2$ echo ${#Zmays_length}
    100
    ljcohen@farm:~/ECL243/Assignment2$ Tdip_length=$(zcat Tdip_R1.fq.gz | sed -n '2p')
    ljcohen@farm:~/ECL243/Assignment2$ echo ${#Tdip_length}
    100


Coverage depth:
* Z mays: 100*1132416/152126075 =  0.7443931
* Z mex: 76*5346216/152126075 = 2.670893
* T dip: 100*17074344/152126075 = 11.22381

###How are the summary stats different?

* Z. mex had 1470419 reads mapped (54.96%) 
* Z. mays had 88127 reads mapped (15.53%)
* T. dip had 602866 read mapped (7.05%)

###Why did so few reads map? What might be an explanaton for any differences?
- Only chromosome 10 from Zea mays was provided as a reference,
- full reference genome might yield better overall coverage
- smallest read length (76bp for Z mex) had highest mapping
- small .fq files provided, more reads might yield better coverage,
- the highest number of reads (Tdip) had the lowest mapping at highest coverage
- species divergence from reference (Tdip)
- Z. mays didn't have very may reads mapping, but only 1.1 million reads input

###samstat
http://rawgit.com/ljcohen/ECL243/master/Assignment2/aln_Tdip.sam.samstat.html
http://rawgit.com/ljcohen/ECL243/master/Assignment2/aln_Zmays.sam.samstat.html
http://rawgit.com/ljcohen/ECL243/master/Assignment2/aln_Zmex.sam.samstat.html

Challenge:
Do you get more or fewer SNPs than you expected? Any idea why?

###samtools reference manual:
http://www.htslib.org/doc/samtools-1.1.html

###bwa reference manual:
http://bio-bwa.sourceforge.net/bwa.shtml

###Helpful links:
https://www.biostars.org/p/638/
http://www.illumina.com/documents/products/technotes/technote_coverage_calculation.pdf
http://arxiv.org/abs/1303.3997
https://www.biostars.org/p/117225/
http://seqanswers.com/forums/showthread.php?t=32708
https://www.biostars.org/p/12475/
https://www.biostars.org/p/84396/
http://davetang.org/wiki/tiki-index.php?page=SAMTools#Simple_stats_using_SAMTools_flagstat
http://davetang.org/wiki/tiki-index.php?page=SAMStat
https://broadinstitute.github.io/picard/explain-flags.html
https://samtools.github.io/hts-specs/SAMv1.pdf

