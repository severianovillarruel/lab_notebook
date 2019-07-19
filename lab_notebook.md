#########################################################Running Velvet: Problem Set 6

/home/svillarr/svillarr/PS/ps6
#7/13/19
1. Wrote the bash script ```num_nt.srun & tot_nt.py``` to calculate the total number of nucleotides that we got back in out reads.
Percent of CPU this job got: 99%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.96
Exit status: 0

800_3_PE5_interleaved.fq_1: 68868840 is nt long
	- paired end
800_3_PE5_interleaved.fq_2: 67942146 is nt long
	- paired end
800_3_PE5_interleaved.fq.unmatched: 62939394 is nt long
	- unpaired end

2. Got number of records:

800_3_PE5_interleaved.fq_1: 858989 records
800_3_PE5_interleaved.fq_2: 858989 records
800_3_PE5_interleaved.fq.unmatched: 849803 records

#7/17/19
/home/svillarr/svillarr/PS/ps6
2. The total genome size for the genome we are assembling is 50 fosmids x 40 Kb long = 2000000
   Meaning the nucleotide coverage for sequencing event is: 68868840/2000000 = 34.43442 coverage at each nucleotide

	 	average contig len = ins_length

    Ck = C * (Lmean - K + 1) / Lmean
    Ck = kmer coverage
    C = coverage = 99.87519
    K = 31, 41, 49
    Lmean = num nt in seq event/num reads in seq event = 199750380/2567781 = 77.7910499377

    Ck = 99.87519 * (77.7910499377 - 31 + 1) / 77.7910499377 = 61.358475
    Ck = 99.87519 * (77.7910499377 - 41 + 1) / 77.7910499377 = 48.51957
    Ck = 99.87519 * (77.7910499377 - 49 + 1) / 77.7910499377 = 38.248446

		```
		Example:
    /usr/bin/time -v velvetg /projects/bgmp/svillarr/PS/ps6/output_dir31/ -exp_cov 61.358475 -ins_length 77.7910499377
		```
		*may not have striped new line characters in ```num_nt.srun & tot_nt.py```. Getting slightly different expected coverage
		length and average contig length*


/home/svillarr/svillarr/PS/ps6
3. Ran velveth with the script I made velveth.srun
	- documentation: http://www.ebi.ac.uk/~zerbino/velvet/
	- need to make an output dir
	*forgot to save standard error*

	kmer 31
		Percent of CPU this job got: 310%
		Elapsed (wall clock) time (h:mm:ss or m:ss): 0:23.47
		Exit status: 0

	kmer 41
		Percent of CPU this job got: 354%
		Elapsed (wall clock) time (h:mm:ss or m:ss): 0:19.82
		Exit status: 0

	kmer 49
		Percent of CPU this job got: 312%
		Elapsed (wall clock) time (h:mm:ss or m:ss): 0:17.79
		Exit status: 0

/home/svillarr/svillarr/PS/ps6
4. Ran on kmers 31, 41, 49 velvetg using velvetg.srun
	 Will run on kmer 49: 20x, 60x, auto, auto_min500nt

  Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir31/ -exp_cov 61.358475 -ins_length 77.7910499377"
	Percent of CPU this job got: 316%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:40.36
	Exit status: 0

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir41/ -exp_cov 48.51957 -ins_length 77.7910499377"
	Percent of CPU this job got: 300%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:31.96
	Exit status: 0

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir49/ -exp_cov 38.248446 -ins_length 77.7910499377"
	Percent of CPU this job got: 280%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:26.54
	Exit status: 0


#7/17/19

	Command exited with non-zero status 1
	Command being timed: "velvetg 0x -cov_cutoff *20x*-ins_length 77.7910499377"
	Percent of CPU this job got: 4%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.04
	Exit status: 1

	velvetg: Could not write to 0x/Log: No such file or directory
	Command exited with non-zero status 1
	Command being timed: "velvetg 0x -cov_cutoff *60x* -ins_length 77.7910499377"
	Percent of CPU this job got: 75%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.00
	Exit status: 1

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir_auto/ -cov_cutoff *auto* -ins_length 77.7910499377"
	Percent of CPU this job got: 265%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:24.25
	Command exited with non-zero status 1

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir_auto_min/ -cov_cutoff *auto -ins_length 77.7910499377 --min_contig_lgth 500"*
	Percent of CPU this job got: 50%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.00
	Exit status: 1

  Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir20x/ -cov_cutoff *20x* -ins_length 77.7910499377"
	Percent of CPU this job got: 265%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:24.26
	Exit status: 0

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir60x/ -cov_cutoff *60x* -ins_length 77.7910499377"
	Percent of CPU this job got: 275%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:22.97
	Exit status: 0

	Command exited with non-zero status 1
	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir_auto_min/ -cov_cutoff *auto -ins_length 77.7910499377 --min_contig_lgth 500"*
	Percent of CPU this job got: 50%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.00
	Exit status: 1

	Command being timed: "velvetg /projects/bgmp/svillarr/PS/ps6/output_dir_auto_min/ -cov_cutoff *auto -ins_length 77.7910499377 -min_contig_lgth 500"*
	Percent of CPU this job got: 255%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:25.22
	Exit status: 0

#7/18/19
	/home/svillarr/svillarr/PS/ps6
	5. Ran ps6_stats.srun & ps6_stats.py on contigs.fa files returned from velvetg
	*converted all the contigs.fa returned from velvetg into one line fasta files using the ```mkfasta_oneline.py  mkfasta_oneline.srun```
	in ~/svillarr/ on talapas*

	kmer 31
		The number of contigs is: 12277
		The total length of the genome assembly is: 2213616 basepairs
		The largest contig is: 3809 basepairs
		The average contig length is: 180.30593793271973 basepairs
		The mean coverage depth is: 36.632877466563684
		The N50 if this genome assembly is: 290


	kmer 41
		The number of contigs is: 5368
		The total length of the genome assembly is: 1753163 basepairs
		The largest contig is: 7443 basepairs
		The average contig length is: 326.59519374068554 basepairs
		The mean coverage depth is: 36.91364448416549
		The N50 if this genome assembly is: 827


	kmer 49
		The number of contigs is: 3172
		The total length of the genome assembly is: 1609227 basepairs
		The largest contig is: 5032 basepairs
		The average contig length is: 507.32250945775536 basepairs
		The mean coverage depth is: 41.30642823770482
		The N50 if this genome assembly is: 1062

	kmer cov 20x
		The number of contigs is: 426
		The total length of the genome assembly is: 396369 basepairs
		The largest contig is: 7268 basepairs
		The average contig length is: 930.443661971831 basepairs
		The mean coverage depth is: 125.69256579577477
		The N50 if this genome assembly is: 1925


	kmer cov 60x
		The number of contigs is: 426
		The total length of the genome assembly is: 396369 basepairs
		The largest contig is: 7268 basepairs
		The average contig length is: 930.443661971831 basepairs
		The mean coverage depth is: 125.69256579577477
		The N50 if this genome assembly is: 1925

	kmer cov auto
		The number of contigs is: 1330
		The total length of the genome assembly is: 1193200 basepairs
		The largest contig is: 10467 basepairs
		The average contig length is: 897.1428571428571 basepairs
		The mean coverage depth is: 58.14230608045102
		The N50 if this genome assembly is: 1796


	kmer cov auto min 500 nt
		The number of contigs is: 620
		The total length of the genome assembly is: 1025835 basepairs
		The largest contig is: 10467 basepairs
		The average contig length is: 1654.5725806451612 basepairs
		The mean coverage depth is: 49.81960170322581
		The N50 if this genome assembly is: 2013

	stats on part 1 fasta
		The number of contigs is: 366
		The total length of the genome assembly is: 864854 basepairs
		The largest contig is: 46013 basepairs
		The average contig length is: 2362.9890710382515 basepairs
		The mean coverage depth is: 392.8630475874315
		The N50 if this genome assembly is: 6819

####################################################Reciprocal Best Hit (RBH) pipeline: Problem Set 7:

#7/11/19
/projects/bgmp/svillarr/PS/ps7
1. Put files on talapas:
  - Homo_sapiens.GRCh38.pep.all.fa
  - Danio_rerio.GRCz11.pep.all.fa
  - human_reference_table.tsv
  - zebrafish_reference_table.tsv

/home/svillarr/svillarr/PS/ps7
2. Converted ``` Homo_sapiens.GRCh38.pep.all.fa & Danio_rerio.GRCz11.pep.all.fa ``` to one line fasta files ```Homo_sapiens.GRCh38.pep.all.oneline.fa & Danio_rerio.GRCz11.pep.all.oneline.fa``` using the ```mkfasta_oneline.py```

3. Wrote a python script to extract the longest proteins associated with Zebrafish and Humans
```
Homo_sapiens.GRCh38.pep.longest.oneline.fa
Danio_rerio.GRCz11.pep.longest.oneline.fa
```

#7/11/19
1. Loaded the module BLAST+/2.2.31 in Talapas
2. Wrote a program called mkfasta_oneline.py to convert multi-line seqs in fasta files to single line seqs


#7/17/19
/home/svillarr/svillarr/PS/ps7/database
1. Made the database.srun script to create a blast database for human and zebrafish longest protein fasta file
2. Made the script palign.srun to align the zebrafish proteins to the human proteins and the human proteins to the zebrafish proteins (should take ~ 7 hours)
 - Slurm Job_id=9615412 Name=protein_align Ended, Run time 10:01:59, COMPLETED, ExitCode 0

local run
./RBH.py -ih qhomo.txt -iz qdanio.txt -th human_reference_table.tsv -tz zebrafish_reference_table.tsv -o output.txt

talapas Run
/usr/bin/time -v ./RBH.py -ih /home/svillarr/svillarr/PS/ps7/database/db_danio.txt -iz /home/svillarr/svillarr/PS/ps7/database/db_homo.txt -th human_reference_table.tsv -tz zebrafish_reference_table.tsv -o output.txt

#########################################################RNA-seq alignment: Problem Set 8: 7/15/19

1. Downloaded: a. Chromosomes: Danio_rerio.GRCz11.dna.chromosome."wildcard".fa.gz
               b. Gene Set :Danio_rerio.GRCz11.97.gtf.gz
   To the talapas directory: /projects/bgmp/svillarr/PS/ps8/dre
2. Installed the short read alignment software STAR version 2.7.1a
3. Installed samtools 1.9
4. Ran the STARprep.srun script on the chromosome and gene set files
        - chrLength.txt, chrName.txt, exonGeTrInfo.tab, geneInfo.tab, genomeParameters.txt, SAindex, sjdbList.fromGTF.out.tab, transcriptInfo.tab, chrNameLength.txt, chrStart.txt, exonInfo.tab, Genome, SA, sjdbInfo.txt, sjdbList.out.tab

      	Percent of CPU this job got: 332%
      	Elapsed (wall clock) time (h:mm:ss or m:ss): 15:13.92
      	Exit status: 0

5. Ran the STARalign.srun script and created the sam file salignAligned.out.sam
        Percent of CPU this job got: 405%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 6:24.74
        Exit status: 0

6. Made the script samtools_script.srun to:
    - view the sam file as a bam file (make a bam file)
    - sort the bam file that was just made
    - index the file that was just sorted
    - view chromosome 1 of the file that was just sorted and put the chromosome 1 data into the file chr1.SAM
