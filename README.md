# Introduction to RNAseq Assembly and Analysis
======
Methods for RNAseq assembly and analysis are not reasonably mature, with several established assemblers and associated analytical pipelines. For a general review of best practices, see Conesa et al. ([2016](http://dx.doi.org/10.1186/s13059-016-0881-8)). Some aspects of existing pipelines are shared by more standard methods for analysis of DNA sequence data. However, working with RNAseq data also presents a number of challenges that may be unfamiliar to researchers who have worked primarily with DNA. One of these challenges is the fact that most genes appear to produce multiple transcript isoforms ([Pan et al. 2008](http://dx.doi.org/10.1038/ng.259)); this presents a challenge to assembly algorithms and results in a far greater diversity of transcripts than one might expect based on estimates of gene diversity. Assemblers and associated applications tend to deal with this challenge by attempting to diagnose clusters of isoforms that are derived from the same gene, but know this will be impossible to do without error given the known difficulty with even distinguishing paralogous copies of the same gene. A second challenge is the fact that we expect RNAseq reads to be present in proportion to their abundance in the tissue from which we obtained RNA; for this reason, most assembly pipelines integrate read frequency data in the QC and analysis processes.

The KU BI Pipeline
======
The pipeline below is designed to be run on the [KU ACF cluster](https://acf.ku.edu/wiki/index.php/Main_Page) and is derived primarily from the recommendations of authors of the most popular de novo assembler for RNAseq data: Trinity. [The Trinity github](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running%20Trinity) includes detailed instructions not only for using the core Trinity assembler, but also for various QC processes associated with this assembler. The main alternative to Trinity involves the use of multiple k-mer assemblers; while Trinity focuses exclusively on 25 k-mer assembly, alternative methods such as trans-ABYSS and trans-SOAP allow users to hierarchically assess the performance of alternative k-mer sizes. While these multiple k-mer approaches can often improve assembly performance, recent comparative studies suggest that the simpler strategy used by Trinity is nearly as effective. Other pipelines available online may also deviate from our recommended Trinity pipeline in one way or another. For example, the [Harvard Pipeline](http://informatics.fas.harvard.edu/best-practices-for-de-novo-transcritome-assembly-with-trinity.html) includes some additional steps for removal of low frequency transcripts or improperly paired reads. In our analyses of whole RNA transcriptomes for an anole lizard, we have generally found that this results in fewer matches with conserved orthologues than we obtain without doing this addtional step and have therefore left this step out of the pipeline below. However, we encourage anyone working with RNAseq to explore these alternatives. Other places to learn more include the [Simple Fools Guide](http://sfg.stanford.edu/quality.html) and the [PASA pipeline](http://pasapipeline.github.io/#A_ComprehensiveTranscriptome).



Step 1: Concatenate Fastq Files
======
Before starting, it may be necessary to concatenate `fastq` formatted sequence output if your transcriptome project was run across different lanes or at different times. This can be accomplished by using the `cat` command to create separate concatenated files for each end of your paired end reads (e.g., one concatenated file for R1 and another for R2). In the example, below, we are concatenating from reads from Illumina runs that were done at two separate times (7June2016 and 19Aug2016) from the same library.
```
cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R1* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R1* >> Digestive_CTTGTA_R1.fastq.gz

cat ../anole_RNAseq_7June2016_run1/Project_Glor_Alexander/Sample_Digestiv/*R2* ../anole_RNAseq_19Aug2016_run2/Project_Glor_Alexander/Sample_Digestiv/*R2* >> Digestive_CTTGTA_R2.fastq.gz

```

Step 2: Preliminary QC & Quality Trimming
======
Preliminary QC of your sequences can be completed by applying [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to your fastq sequence files. [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a popular tool that "aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines." Carefully inspect the `html` output from [`fastqc`](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) because this is the main way that you are going to make sure that nothing is seriously wrong with your data before you delve into the time consuming series of analyses discussed below.
```
fastqc -k 6 *
```

If your sequences look OK after preliminary QC, its time to get your sequences ready for downstream analyses by trimming adaptors and eliminating low quality sequences and basecalls. We are going to do this by using the function `cutadapt` to (1) trim Illumina adapter sequences, (2) discard reads <25 bp in length (otherwise *de novo* assembly in Trinity will fail because its kmer = 25) and (3) perform gentle trimming of low quality basecalls. Recent studies suggest that trimming to a relatively low PHRED score of 5 results in transcriptomes that are considerably more complete than those that result from more aggressive quality trimming, without a commensurate increase in errors ([MacManes 2014](http://journal.frontiersin.org/article/10.3389/fgene.2014.00013/full)). Running cut adapt for datasets the size of most of our individual transcriptomes is computationally non-trivial and the function does not use multiple threads, so best to give run your job with a walltime of 6 hours for datasets with 30-50 million reads.
```
#PBS -N cutadapt.sh
#PBS -l nodes=1:ppn=1:avx,mem=16000m,walltime=5:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Testes
#PBS -j oe
#PBS -o cutadapterror

cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 5 -m 25 -o Testes_trimmed_R1.fastq.gz -p Testes_trimmed_R2.fastq.gz Testes_CTTGTA_R1.fastq.gz Testes_CTTGTA_R2.fastq.gz > cutadapt.log
```
After trimming is complete, use `fastqc` on the resulting files to check that adapters have been trimmed and that the newly generated `fastq` files look good. 

Step 3: Assembly
======
After you've conducted the basic QC steps described above, you're ready to do assemble your transcriptome. We are going to assemble transcriptomes using three methods: (1) *de novo* assembly, (2) genome guided assembly, and (3) transcriptome guided assembly.

Step 3a: *de novo* Assembly with Trinity
------
After you've conducted the basic QC steps described above, you're ready to do your first *de novo* assembly. We use the `Trinity` package for *de novo* transcriptome assembly. Because *de novo* assembly with `Trinity` requires a large amount of computer memory (~1GB RAM/1 million sequence reads) and generates a large number of files (~900,000) you should be sure that your quotas are able to accommodate these requirements prior to initiating assembly. In the script below, we will run *de novo* assembly using a high memory node via the `bigm` queue; moreover, hundreds of thousands of files that are temporarily required during the assembly process are stored on the node running the analyses (`file=200gb`) rather than in the scratch space, which has stricter constraints on file size and number. Importantly, all of the temporary files generated during the course of the `Trinity` analyses are deleted at the end of the script to avoid gumming up the cluster's hard drives. The main output from `Trinity` will be a large (e.g., 400+ MB for a dataset with 40 million reads) `fasta` file called `trinity_out_dir.Trinity.fasta` that will include each of your assembled transcripts. This operation will take a day or more running on multiple threads.

```
#PBS -N trinity_heart
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=512000m,walltime=72:00:00,file=200gb
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart
#PBS -j oe
#PBS -o trinity_heart_error

work_dir=$(mktemp -d) #Generates name for working directory where temporary files will be stored.
mkdir $work_dir #Generates working directory for temporary file storage.
cp /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart/Heart_trimmed* $work_dir #Copies sequence files for assembly to new working directory on node doing analyses.
Trinity -seqType fq --max_memory 512G --left $work_dir/Heart_trimmed_R1.fastq.gz --right $work_dir/Heart_trimmed_R2.fastq.gz -SS_lib_type RF --CPU 24 --full_cleanup --normalize_reads --min_kmer_cov 2 --output $work_dir/trinity_out_dir >> trinity.log
mv $work_dir/trinity_out_dir.Trinity.fasta /scratch/glor_lab/rich/distichus_genome_RNAseq/Heart/ #Moves files from node doing analyses to scratch drive.
rm -rf $work_dir

```

Step 3b: Reference Genome-based Assembly via TopHat and Cufflinks
-----
An alternative to *de novo* assembly is to assemble transcriptomes using a previously assembled genome as a reference. This procedure can be implemented via the Tuxedo software suite, wherein the program Tophat is first used to map the individual Illumina sequence reads to the genome and to identify putative splice sites. Once this mapping is complete, Cufflinks may be used to assemble transcripts against the genome.

Tophat uses the mapping function Bowtie, but differs in being able to handle the types of larger gaps that are expected in transcriptomes where introns and alternative splicing may result in gaps of hundreds or even thousands of basepairs. The output from Tophat will be in the form of a large BAM file that includes the mapping information. Tophat I'm not sure why, but the Tophat function below was incredibly time and storage intensive, ultimately requiring more than two weeks to complete and requiring more than 400 GB of temporary storage.

```
#PBS -N tophat_dovetail_short
#PBS -q bigm -l nodes=1:ppn=24:avx,mem=500000m,walltime=296:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Cufflinks
#PBS -j oe
#PBS -o tophat_dovetail_short_error

bowtie2-build /scratch/a499a400/anolis/dovetail/scrubbed_genome.fasta Dovetail_bowtie_DB.fasta
tophat -r 0 Dovetail_bowtie_DB.fasta /scratch/glor_lab/rich/distichus_genome_RNAseq/All_Tissues/decontam_RNA_1.fastq.gz /scratch/glor_lab/rich/distichus_genome_RNAseq/All_Tissues/decontam_RNA_2.fastq.gz
```
Once we have the desired BAM file from Tophat, we can use Cufflinks to generate transcripts as follows:
```
#PBS -N cufflinks_all
#PBS -q default -l nodes=1:ppn=24:avx,mem=50000m,walltime=148:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome/Cufflinks
#PBS -j oe
#PBS -o cufflinks_all_error

cufflinks /scratch/glor_lab/rich/distichus_genome/Cufflinks/tophat_out/accepted_hits.bam
```
Step 3b: Reference Transcriptome-based Assembly via BWA
-----
Here we will use BWA to assemble our transcriptome with reference to the previously published and curated transriptome for *Anolis carolinensis*. 
```
BWA_ALN1="bwa aln -n 5 -q 20 -t 4 -f $SAI_FILE1 $REF_FILE /home/data/anoles/synced_clt_${NAMES[$i]}_R1.fastq"
        BWA_ALN2="bwa aln -n 5 -q 20 -t 4 -f $SAI_FILE2 $REF_FILE /home/data/anoles/synced_clt_${NAMES[$i]}_R2.fastq"
BWA_SAMPE="bwa sampe -P $REF_FILE $SAI_FILE1 $SAI_FILE2 /home/data/anoles/synced_clt_${NAMES[$i]}_R1.fastq /home/data/anoles/synced_clt_${NAMES[$i]}_R2.fastq 2> /$

        SAMTOOLS_SORT="samtools sort $UNSORT_BAM_FILE ${NAMES[$i]}.pre"

        SAMTOOLS_CALMD="samtools calmd -AEbr ${NAMES[$i]}.pre.bam $REF_FILE 2> /dev/null > ${NAMES[$i]}.bam"
```

Step 4: Basic Transcriptome QC
======
We will conduct several QC analyses with the Trinity assembly, most of which are directly recommended by the authors of Trinity. In this step, we will [map our reads back to the assembly](https://github.com/trinityrnaseq/trinityrnaseq/wiki/RNA-Seq-Read-Representation-by-Trinity-Assembly) with the expectation that most (>70%) will be successfully mapped. We will then [calculate our transriptome N50](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome%20Contig%20Nx%20and%20ExN50%20stats). 

Step 4a: Mapping Reads to Assembled Transcriptome
------
Our first assembly QC step will involve mapping our individual sequencing reads back to the assembly. This step is necessary to ensure that most reads (70-80%) can be mapped to the assembly as "proper pairs;" those that do not are likely the result of extremely low frequency, erroneous sequences, or instances where one of the two paired end sequences failed. In order to complete this operation, we must have access to [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) and [samtools](http://samtools.sourceforge.net/), both of which are available via the `env-selector-menu` option on the cluster. Once all the appropriate software is available, we need to build a `bowtie2` index. This operation will produce six files with the extension `bt2` and should take 10-15 minutes run on a single processor in interactive mode.
```
bowtie2-build trinity_out_dir.Trinity.fasta Trinity_testes.fasta
```
Next we do the actual mapping, which is best submitted to the default queue because this function may take a day or more to complete. I don't know how parallelized this process is but assign it to 24 processors.
```
#PBS -N mapping_testes
#PBS -q default -l nodes=1:ppn=24:avx,mem=50000m,walltime=24:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Testes
#PBS -j oe
#PBS -o mapping_testes_error
bowtie2 --local --no-unal -x trinity.fasta -q -1 Testes_trimmed_R1.fastq.gz -2 Testes_trimmed_R2.fastq.gz | samtools view -Sb - | samtools sort -no - - > bowtie2.nameSorted.bam
```
To get a summary of the mapping results, we can use a perl script from `Trinity`.
```
SAM_nameSorted_to_uniq_count_stats.pl bowtie2.nameSorted.bam

#read_type      count   pct
proper_pairs    30446172        94.39
improper_pairs  928085  2.88
left_only       699683  2.17
right_only      181046  0.56

Total aligned rnaseq fragments: 32254986
```

Step 4b: Transcriptome Nx Statistic
------
In this step, we calculate basic statistics to assess contig assembly. Nx is the transcriptome equivalent of N50. 

```
TrinityStats.pl trinity_out_dir.Trinity.fasta > transcriptome_contig_nx_stat.log
```
The authors of Trinity encourage assessment of Nx based on the single longest transcript per isoform given that the tendency of assemblers to produce too many isoforms, particularly for larger transcripts (as a result, we paradoxically expect that our N50 for the dataset with only the longest isoform will actually be lower than that for all transcripts).

```
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes':  456614
Total trinity transcripts:      590779
Percent GC: 42.48

########################################
Stats based on ALL transcript contigs:
########################################

        Contig N10: 3701
        Contig N20: 2409
        Contig N30: 1623
        Contig N40: 1072
        Contig N50: 703

        Median contig length: 304
        Average contig: 535.45
        Total assembled bases: 316335067


#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

        Contig N10: 3065
        Contig N20: 1523
        Contig N30: 837
        Contig N40: 561
        Contig N50: 426

        Median contig length: 283
        Average contig: 425.57
        Total assembled bases: 194318977
```

Step 5: Preliminary Estimates of Transcript Abundance and Transcriptome Ex
======
N50 is a very important statistic in genomics, but does not capture all the information we might like to consider when assembling transcriptomes. For this reason, transcriptome QC also involves calculation of Ex, a variant of Nx that incoporporations information on transcript frequency. In order to calculate Ex, we must first estimate transcript abundance using 'RSEM' ([Teng et al. 2016](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4842274/). We will use `RSEM` via a script in the `Trinity` toolkit, but we still need to have `RSEM` in our path before proceeding. Before running the functions below, it is also necessary to have access to Bowtie and express. The following command will prepare the reference, generate an alignment and estimate transcript abundance in interactive mode; it should only take 10 or 15 minutes to run, but will generate at least 1GB of files with `trinity_out_dir.Trinity.fasta.bowtie` or `trinity_out_dir.Trinity.fasta.RSEM` prefixes.

```
align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left Testes_trimmed_R1.fastq.gz --right Testes_trimmed_R2.fastq.gz --SS_lib_type RF --thread_count 24 --est_method RSEM --output_dir trin_rsem --aln_method bowtie --trinity_mode --prep_reference
```
We should now be able to construct matrices containing counts and a matrix of normalized expression values for isoforms (with `trans_counts` prefix)  and genes (with `gene_counts` prefix). Normalized values will be calculated as transcripts per million transcripts (TPM) and FPKM. Running these perl scripts included in the `Trinity` toolkit should only take a few seconds in interactive mode.
```
abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix trans_counts --name_sample_by_basedir trin_rsem/RSEM.isoforms.results
abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix gene_counts --name_sample_by_basedir trin_rsem/RSEM.genes.results
```
The resulting matrix containing normalized expression values won't take up too much space.
```
transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
TRINITY_DN0_c0_g1_i1    TRINITY_DN0_c0_g1       224     64.93   2.00    0.80    1.08    100.00
TRINITY_DN100000_c0_g1_i1       TRINITY_DN100000_c0_g1  237     76.47   5.71    1.93    2.62    100.00
TRINITY_DN100000_c0_g2_i1       TRINITY_DN100000_c0_g2  237     76.47   2.29    0.77    1.05    100.00
TRINITY_DN100001_c0_g1_i1       TRINITY_DN100001_c0_g1  265     102.34  0.00    0.00    0.00    0.00
TRINITY_DN100002_c0_g1_i1       TRINITY_DN100002_c0_g1  367     201.70  4.50    0.58    0.78    100.00
TRINITY_DN100002_c0_g2_i1       TRINITY_DN100002_c0_g2  367     201.70  4.50    0.58    0.78    100.00
TRINITY_DN100004_c0_g1_i1       TRINITY_DN100004_c0_g1  269     106.12  3.00    0.73    0.99    100.00
TRINITY_DN100006_c0_g1_i1       TRINITY_DN100006_c0_g1  339     174.02  10.60   1.57    2.14    76.47
TRINITY_DN100006_c0_g1_i2       TRINITY_DN100006_c0_g1  292     128.12  2.40    0.48    0.66    23.53
```

We can now obtain a summary of how many genes are expressed at different normalized expression values. This is useful for determining how many transcripts are present at different frequencies and potentially determining how many of your transcripts are very lowly expressed and possibly erroneous. Note that the required function is part of the `Trinity` toolkit, but the folder that this function is in is not included as part of your path when you add `Trinity` to your work environment on the cluster, which explains why I have used explicit directory infomration below (if you don't do this or ad the `misc` folder in the `Trinity` `util` folder to your path you will get an error).
```
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl trans_counts.TPM.not_cross_norm > trans_matrix.TPM.not_cross_norm.counts_by_min_TPM

/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl gene_counts.TPM.not_cross_norm > genes_matrix.TPM.not_cross_norm.counts_by_min_TPM
```
The resulting matrix of numbers should look something like the table below, where the first column indicates normalized expression and the second column indicates the number of transcripts with that level of expression. In the table below, for exampe, 11,904 transcripts are expessed with a TPM of 5 whereas 341953 are expressed with a TPM of 1. 

```
neg_min_tpm     num_features
-108409 1
-51764  2
-35403  3
-32386  4
-21368  5
...   ...
-5      11904
-4      16134
-3      26908
-2      78673
-1      284410
0       341953
```

We can use the table above to generate a plot of expression versus numbers of transcripts; such plots can be useful for figuring out how many genes are well-supported by your expression data. Instructions for going this process in R are available via `Trinity` and are replicated below.

```
data = read.table("Dewlap_genes_matrix.TPM.not_cross_norm.counts_by_min_TPM", header=T)
plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')
filt_data = data[data[,1] > -100 & data[,1] < -10,]
fit = lm(filt_data[,2] ~ filt_data[,1])
print(fit)
abline(fit, col='green', lwd=3)
```

The y-intercept of this fitted line, which is ~2,500 in the case of the sample used here, can be interpreted as a somewhat more reasonable estimate of the genes expressed in a sample than prior estimates based on the number of assembled transcripts.

```
Call:
lm(formula = filt_data[, 2] ~ filt_data[, 1])

Coefficients:
   (Intercept)  filt_data[, 1]  
       2493.53           26.55
```
In addition to using the transcript abundance data to generate this rough estimate of how many genes are well-supported by your data, we will also generate the Ex90N50 value which is calculated in the same manner as N50, but only includes samples above some threshold of transcript frequency. For example, E90 includes only those samples included in the 90% transcript frequency. To generate Ex90N50 estimates, go into `trin_rsem folder`, and run the following commands to extract the columns of interest (because we are only working on one sample)
```
cat RSEM.isoforms.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.isoforms.results.mini_matrix
cat RSEM.genes.results  | perl -lane 'print "$F[0]\t$F[5]";' >  RSEM.genes.results.mini_matrix
```
Then move back to the folder containing the trinity.fasta file and run the following commands.
```
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.isoforms.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_trans.stats
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/contig_ExN50_statistic.pl trin_rsem/RSEM.genes.results.mini_matrix trinity_out_dir.Trinity.fasta > ExN50_genes.stats
```
The resulting table (below) indicates that 90% of the transcripts are accounted for by 218,896 transcripts with an N50 of 896.

```
#E      min_expr        E-N50   num_transcripts
E10     108408.66       354     1
E16     51763.51        1902    2
E19     35391.32        1658    3
E21     24329.59        1902    4
E24     21351.73        1902    5
E25     16020.57        1902    6
E27     15276.74        1902    7
...
E90     0.89    896     218896
E91     0.85    869     230431
E92     0.81    846     242485
E93     0.77    822     255135
E94     0.73    799     268501
E95     0.68    778     282737
E96     0.62    762     298121
E97     0.56    745     315158
E98     0.49    730     334326
E99     0.39    721     356997
E100    0       742     439159
```

Step 6: Assess Full-length Transcript Coverage 
======
One important benchmark for transcriptome quality involves asking how many assembled transcripts correspond with full-length transcripts in other organisms. Such counts can be generated in several ways, most frequently using `blastx` and `BUSCO` to compare transcripts to a database of the users choosing or a curated database of conserved orthologues, respectively. This operation will take at 20+ hours to complete.

Step 6a: Benchmarking Universal Single-Copy Orthologs (BUSCO)
------
"BUSCO v2 provides quantitative measures for the assessment of genome assembly, gene set, and transcriptome completeness, based on evolutionarily-informed expectations of gene content from near-universal single-copy orthologs selected from OrthoDB v9." For more about BUSCO visit the [project's website](http://busco.ezlab.org/) or read the paper reporting the method ([Simão et al. 2015](http://dx.doi.org/10.1093/bioinformatics/btv351)).
```
#PBS -N busco_dewlap.sh
#PBS -l nodes=1:ppn=24:avx,mem=64000m,walltime=72:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Dewlap
#PBS -j oe
#PBS -o busco_dewlap_error

BUSCO.py -o busco_dewlap -i trinity_out_dir.Trinity.fasta -l ../vertebrata_odb9 -m tran
```
After running this operation, you should find a simple summary of your data (`short_summary_busco_dewlap.txt`) in a folder called `run_busco_name` that will tell you how many complete BUSCOs were recovered in your dataset.
```
        1804    Complete BUSCOs (C)
        1006    Complete and single-copy BUSCOs (S)
        798     Complete and duplicated BUSCOs (D)
        435     Fragmented BUSCOs (F)
        347     Missing BUSCOs (M)
        2586    Total BUSCO groups searched
```

Step 6b: Assess Full-length Transcripts Relative to Reference Via BLAST+
------
A second approach to determining how many assembled transcripts correspond with transcripts in previously sequenced datasets is to use BLAST+ to queuery your transcripts against a database of previously identified transcripts. If you do not have a reference transcriptome or proteome from your focal species (or anything close to your species), you can conduct this comparative analysis using a more general reference from SwissProt or elsewhere. The Trinity GitHub [provides details on conducting this analysis](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts), which are revised below for use at the KU ACF.

```
makeblastdb -in ../uniprot_sprot.fasta -dbtype prot
blastx -query trinity_out_dir.Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6
analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fast > analyze_blastPlus_topHit_coverage.log
blast_outfmt6_group_segments.pl blastx.outfmt6 trinity_out_dir.Trinity.fasta uniprot_sprot.fasta > blast.outfmt6.grouped
blast_outfmt6_group_segments.tophit_coverage.pl blast.outfmt6.grouped > analyze_groupsegments_topHit_coverage.log

```
With anoles, we can use the proteome generated from the Anolis carolinensis genome project. This should take less than 24 hours.

```
#PBS -N blast_anole_skin.sh
#PBS -l nodes=1:ppn=12:avx,mem=200000m,walltime=24:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Skin/Anole_BLAST
#PBS -j oe
#PBS -o blast_anole_skin_error

makeblastdb -in ASU_Acar_v2.1.prot.fa -dbtype prot
blastx -query ../trinity_out_dir.Trinity.fasta -db ASU_Acar_v2.1.prot.fa -out blastx.outfmt6 -evalue 1e-20 -num_threads 12 -max_target_seqs 1 -outfmt 6
analyze_blastPlus_topHit_coverage.pl blastx.outfmt6 ../trinity_out_dir.Trinity.fasta ASU_Acar_v2.1.prot.fa > analyze_blastPlus_topHit_coverage.log
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/blast_outfmt6_group_segments.pl blastx.outfmt6 ../trinity_out_dir.Trinity.fasta ASU_Acar_v2.1.prot.fa > blast.outfmt6.grouped
/tools/cluster/6.2/trinityrnaseq/2.2.0/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blast.outfmt6.grouped > analyze_groupsegments_topHit_coverage.log
```

Step 7: Generate TransRate Score
======
TransRate is a program for assessing overall transcriptome quality. Running Transrate should only take about 30 minutes to more than an hour if run across 12 processors.
```
#PBS -N transrate_eye.sh
#PBS -l nodes=1:ppn=12:avx,mem=200000m,walltime=2:00:00
#PBS -M glor@ku.edu
#PBS -m abe
#PBS -d /scratch/glor_lab/rich/distichus_genome_RNAseq/Eye
#PBS -j oe
#PBS -o transrate_eye_error

transrate --assembly trinity_out_dir.Trinity.fasta --left Eye_trimmed_R1.fastq.gz --right Eye_trimmed_R2.fastq.gz --threads 1 >> transrate_eye.log
```
TransRate (as of the 15-Oct-2016), transrate has some issues as mentioned in https://github.com/blahah/transrate/issues/201, so running both v1.0.1 and v.1.0.3 and splicing together their output should be attempted. 
```
`transrate --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate103.log

mv transrate_results transrate_results_103

transrate _1.0.1_ --assembly trinity_out_dir.Trinity.fasta --left Brain_trimmed_R1.fastq.gz --right Brain_trimmed_R2.fastq.gz --threads 8 >> transrate101.log
```

Step 8: Compute DETONATE scores
======
```
#RSEM (reference-free mode)
/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-estimate-transcript-length-distribution trinity_out_dir.Trinity.fasta /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt

/public/detonate-1.11-precompiled/rsem-eval/rsem-eval-calculate-score --paired-end --bam trin_rsem/bowtie.bam trinity_out_dir.Trinity.fasta sample_brain 200 --transcript-length-parameters /public/detonate-1.11-precompiled/rsem-eval/true_transcript_length_distribution/anolis_distichus.txt --strand-specific -p 4  >& rsem_eval.log
```



