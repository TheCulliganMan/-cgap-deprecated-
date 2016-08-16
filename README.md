# cgap
##What it does:
Extracts genes using a consensus sequence.
##How it works:
1. Creates local short read databases from fastq file for blast.
2. Blasts reference fastas against the short read blast databases.
3. Collects blast hits into .hits files.
4. Creates new small fastq files from .hits files.
5. Uses bwa to create gene alignments.
6. Uses bcftools, samtools, and novocraft to create consensus sequences.
##Note
This is still very much a work in progress.  The cgap file needs to be split into multiple functions instead of the monolithic class in which it currently resides.
