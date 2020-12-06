# Viral Sequencing Reads Analysis Piplneline

### Xinpei Zhang

## CONTENT
- Preperation
  - Remote Server Instruction
  - Configurate conda environment
  - Software installation
- Procedure
  - miniconda
  - fastp
  - bwa-mem
  - smtools
  - sambamba
  - SPAdes
  - metaquast
  - pilon
  - mafft
  - iqtree
  - iTOL
- 3.Questions



## Prepreration

### Remote server Instruction

#### 1.Download the remote connection software

Download **Xshell** from [NetSarang.com][1] through `>Free for family/school` manner.

![xshell.png](https://github.com/Simple53/Markdown_Photos/blob/main/xshell.png?raw=true) 

#### 2.Set up a new account by fill in the authentical information.

![ssh.png](https://github.com/Simple53/Markdown_Photos/blob/main/ssh.png?raw=true)



#### 3.Node Changing

`ssh node1 #login to the computing node` 

**Notice!**: `tc6000` is a networked login node, `node1` is a computing node that is unnetworked. Thus, we cannot run a job on tc6000 node, otherwise the server will crash. 

###Configurate conda environment 

Miniconda is a free minimal installer for conda. It is a small, bootstrap version of Anaconda that includes only conda, Python, the packages they depend on, and a small number of other useful packages, including pip, zlib and a few others. Use the `conda install` to install 720+ additional conda packages from the Anaconda repository. 

#### - Install miniconda

**Download miniconda**

 ```bash
wget -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-latest-Linux-x86_64.sh # download through Tsinghua channel
bash Miniconda2-latest-Linux-x86_64.sh # run installation package
# Type Enter or yes
 ```

**Check the install condition**

```
source ~/.bashrc #Re-execute the latest modified initialization file
conda --help 
```
**Add download source**
```bash
conda config --add channels r
conda config --add channels conda-forge
conda config --add channels bioconda
cat  ~/.condarc
```
**Create environment**

```bash
conda create -n env_name python=2
```
**Basic Instructions**

- Entry environment
  `source activate env_name` or `conda activate env_name`
- Exit environment
  `deactivate env_name `
- Delete envirnment
  `conda remove -n env_name --all` or
  `conda env remove  -n env_name`
- View all the environments 
`conda info -e`
- Upgrade all packages
`conda upgrade pack_name`
- Install packages
`conda install pack_name`
- View all 
`conda list`
- Delete packages
`conda remove pack_name`

#### - Software Installation

**Package List**

> fastp，fastqc，trim-galore，bwa，sambamba，samtools，SPAdes，quast，pilon，mafft，qtree

**Install command**

```bash
  conda install [-c bionconda] -y pack_name  
  # -c Additional channel to search for packages
  # -y Do not ask for confirmation
```

## Procedure

#### 1. Fastp——cut adapters and filter out bad reads
```bash
fastp -i ${Indir}/seq_1.fq.gz\
-I ${Indir}/seq_2.fq.gz\
-o ${Outdir}/seq_1.fq.gz\
-O ${Outdir}/seq_2.fq.gz  
```
#### 2. Bwa mem——reads alignment

BWA（Burrows-Wheeler Alignment） is a fast light-weighted tool that aligns short sequences to a sequence database, such as the human reference genome. By default, BWA finds an alignment within edit distance 2 to the query sequence, except for disallowing gaps close to the end of the query. It can also be tuned to find a fraction of longer gaps at the cost of speed and of more false alignments. 

- **Build index**

  `bwa index $Indir/ref.fasta -p $Outdir/ref -a is`

  | **Arguments** | Contents                          |
  | ------------- | --------------------------------- |
  | -p            | prefix of the index               |
  | -a            | Algorithm model: bwtsw, is or rb2 |

- **Align reads**

```
bwa mem -M -t 10 -R "@RG\tID:\tLB:HMvero\tPL:ILLUMINA\tSM:Cov" \
Seq_1.fq.gz  Seq_2.fq.gz > $Outdir/seq.sam 
```

| Arguments | Contents                                             |
| :-------- | ---------------------------------------------------- |
| -M        | mark shorter split hits as secondary                 |
| -t        | number of threads                                    |
| -R        | read group header line such as '@RG\tID:foo\tSM:bar' |
#### 3. Samtools

**3.1 Samtools view——filter and convert .sam to .bam**

```bash 
samtools view -@ 10 -F 12 -bhS \
$Indir/align.sam > $Outdir/align.bam
```

| Arguments | Contents                                                     |
| :-------- | ------------------------------------------------------------ |
| -@        | number of threads                                            |
| -F        | Filter options——check on [Picard](https://broadinstitute.github.io/picard/explain-flags.html)<br />The flag value will be shown in the SAM Flag field. |
| -b        | output BAM                                                   |
| -h        | include header in SAM output                                 |
| -S        | ignored (input format is auto-detected)                      |


**3.2 Samtools flagstat——counts the number of alignments for each FLAG type** 

```bash
samtools flagstat -@ 4 \
$Outdir/align.sam >  $Outdir/align.txt
```

Output Example:

```shell
479111218 + 0 in total (QC-passed reads + QC-failed reads)
3439348 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
473770860 + 0 mapped (98.89% : N/A)
475671870 + 0 paired in sequencing
237835935 + 0 read1
237835935 + 0 read2
462393464 + 0 properly paired (97.21% : N/A)
469984246 + 0 with itself and mate mapped
347266 + 0 singletons (0.07% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Propotion of the mapped reads(Sample_1): 98.89%

Explanation of each line in result file：

1. The number of QC passed reads, QC failed reads, Total reads;
2. The number of repeated reads, QC pass and Failed;
3. The number of reads compared to the reference genome;
4. Quantity of paired reads data;
5. The number of READ1;
6. The number of READ2;
7. The number of reads that correctly match the reference sequence;
8. A pair of reads were aligned to the number of reference sequences, but they were not necessarily aligned to the same chromosome;
9. In a pair of reads, there is only one number matching the reference sequence; 
10. The number of different chromosomes was detected with a pair of reads;
11. A pair of reads were aligned to the number of different chromosomes with a mass value greater than 5.

**3.3 Smatools sort——sorts SAM/BAM/CRAM files**

Sort alignments by leftmost coordinates, or by read name when **-n** is used. 

```bash
samtools sort -@ 20 $Indir/In.bam -o $Outdir/out_sort.bam 
```

| Arguments | Contents                                                     |
| :-------- | ------------------------------------------------------------ |
| -@        | Number of threads                                            |
| -b        | Create a BAI index. This is currently the default when no format options are used. |
| -c        | Create a CSI index. By default, the minimum interval size for the index is 2^14, which is the same as the fixed value used by the BAI format. |
| -m        | Create a CSI index, with a minimum interval size of 2^INT    |

**3.4 Samtools index——Index for following analysis** 

```bash
samtools index -b -@ 20 $Indir/sort.bam 
```

#### 4. Sambamba——subsample bam files

Sambamba allows to efficiently filter BAM file for alignments satisfying various conditions, as well as access its SAM header and information about reference sequences. In order to make these data readily available for consumption by scripts in Perl/Python/Ruby, JSON output is provided. 

- Caculate the coverage：

$$
Coverage=reads * length / genomesize)
$$

- Subsample your reads below 10,000 X coverage
  $$
  parameter≈10000/coverage
  $$




```bash
sambamba view -t 20 -f bam -s 0.01 \
$Indir/sort.bam -o $Outdir/subsample.bam
```

| Arguments | Contents              |
| :-------- | --------------------- |
| -t        | Number of threads     |
| -f        | File format           |
| -s        | Fraction of subsample |

- Use samtools to convert bam to fastq file

  ```bash
  samtools fastq -@ 20 \
  $Indir/subsample.bam \
  -1 $Outdir/subsample_1.fq.gz \
  -2 $Outdir/subsample_2.fq.gz \
  -0 /dev/null -n 
  ```

| Arguments | Contents                               |
| :-------- | -------------------------------------- |
| -@        | Number of threads                      |
| -0        | Dumped reads                           |
| -1        | Write to read1 of paired reads         |
| -2        | Write to read2 of paired reads         |
| -n        | Don't append /1 or /2 to the read name |

#### 5.SPAdes——genome assembly

SPAdes is a program developed by The St. Petersburg Academic University of the Russian Academy of Sciences in collaboration with U.S. scientists to splice data from small genomes, such as bacteria and fungi.The latest version supports the common Illumina Miseq/Hiseq and Ion Torrent sequencing data, and the sequencing data of pacBio and nanopore on the single-molecule sequencing platform can also be assembled, as well as the mixed data.The GAGE-B results on the Miseq platform received the best evaluation. 

```bash
spades.py 
-t 20 --[model]\
-k 21,33,55,77,99,127 \
--12 $Indir/subsample.fastq \
-o $Outdir/spades 
```

| Arguments | Contents                                                     |
| :-------- | ------------------------------------------------------------ |
| -t        | number of threads                                            |
| -model    | --sc this flag is required for MDA (single-cell) data<br />--meta this flag is required for metagenomic sample data<br />--isolate  this flag is highly recommended for high-coverage isolate and multi-cell data |
| -k        | list of k-mer sizes (must be odd and less than 128)          |
| -12       | file with interlaced forward and reverse paired-end reads    |
| -o        | directory to store all the resulting files (required)        |

#### 6.Meta-quast——evaluate the assemble quality 

MetaQUAST evaluates and compares metagenome assemblies based on alignments to close references. It is based on [QUAST](http://cab.spbu.ru/software/quast/) genome quality assessment tool, but addresses features specific for metagenome datasets:

- *Huge species diversity* – the tool accepts multiple references and makes multi-genome tables and plots, including Krona charts.
- *Commonly unknown species content* – the tool automatically detects and downloads reference sequences from NCBI.
- *Presence of highly relative genomes* – the tool detects chimeric contigs and reports “interspecies misassemblies” in addition to the regular assembly error types.

```bash 
metaquast 
-t 20 \
-r  ref.fasta \
test1.fasta test2.fasta [test3,test4,……]\
-o $Outdir/compare
```

| Arguments | Contents                                              |
| :-------- | ----------------------------------------------------- |
| -t        | number of threads                                     |
| -r        | reference fasta                                       |
| -o        | directory to store all the resulting files (required) |
|           |                                                       |

**Results**

From the statitic data we could conclude that the `--meta` approach performs the best. Thus ,  the following  analysis will be carried  on using the meta  fasta file.

![quast.png](https://github.com/Simple53/Markdown_Photos/blob/main/quast.png?raw=true) 



#### 7. Pilon

![皮隆](https://camo.githubusercontent.com/af1086f76a563982f2df133b9a168026de7daa4a1d91c75cfdfc748a718741ff/687474703a2f2f73332e616d617a6f6e6177732e636f6d2f70696c6f6e2f70696c6f6e5f6c6f676f2e706e67)

Pilon is a software tool which can be used to Automatically improve draft assemblies and find variation among strains. Pilon requires as input a FASTA file of the genome along with one or more BAM files of reads aligned to the input FASTA file. Pilon uses read alignment analysis to identify inconsistencies between the input genome and the evidence in the reads. It then attempts to make improvements to the input genome, including:

- Single base differences
- Small indels
- Larger indel or block substitution events
- Gap filling
- Identification of local misassemblies, including optional opening of new gaps

Pilon then outputs a FASTA file containing an improved representation of the genome from the read data and an optional VCF file detailing variation seen between the read data and the input genome.

Analysis Process:

Step1：Alignment

```
bwa index -p index/draft draft.fa
bwa mem -t 20 index/draft read_1.fq.gz read_2.fq.gz | samtools sort -@ 10 -O bam -o align.bam
samtools index -@ 10 align.bam
```

Step2：Mark duplicates

```
sambamba markdup -t 10 align.bam align_markdup.bam
```

Step3：Filter reads

```
samtools view -@ 10 -q 30 align_markdup.bam > align_filter.bam
samtools index -@ 10 align_filter.bam
```

Step4：Polish through Pilon

```
MEMORY= # Depend on the genome size
java -Xmx${MEMORY}G -jar pilon-1.23.jar \
--genome draft.fa --frags align_filer.bam \
--fix snps,indels \
--output pilon_polished --vcf 
```

Notice：The pilon treatment needs to be repeat for an ideal consequence,but remember doing no more than three times.

| Arguments | Contents                                                     |
| :-------- | ------------------------------------------------------------ |
| -Xmx      | Setup enough memory for pilon                                |
| -jar      | run java files                                               |
| --genome  | reference fasta                                              |
| --frags   | A bam file consisting of fragment paired-end alignments      |
| --fix     | A comma-separated list of categories of issues to try to fix |
| --output  | prefix of the resulting files                                |
| --outdir  | dictionary of the final files                                |


#### 8. Mafft

MAFFT is a program used to create multiple sequence alignments of amino acid or nucleotide sequences. Published in 2002, the first version of MAFFT used an algorithm based on progressive alignment, in which the sequences were clustered with the help of the Fast Fourier Transform. Subsequent versions of MAFFT have added other algorithms and modes of operation, including options for faster alignment of large numbers of sequences, higher accuracy alignments, alignment of ncRNA sequences, and the addition of new sequences to existing alignments.

```bash
mafft --thread 30 \
--reorder $Indir/pilon.fasta > Outdir/mafft.fasta
```

| Arguments | Contents                                |
| :-------- | --------------------------------------- |
| --thread  | Number of threads                       |
| --reorder | Outorder: aligned, default: input order |

#### 9. IQ-tree

IQ-TREE software was created as the successor of IQPNNI and TREE-PUZZLE. IQ-TREE can handle a large amount of data and provide more complex models of sequence evolution. As output IQ-TREE will write a self-readable report file `.iqtree`, a NEWICK tree file `.treefile` , which can be visualized by viewer programs such as [FigTree](http://tree.bio.ed.ac.uk/software/figtree/), [Dendroscope](http://dendroscope.org/) or [iTOL](http://itol.embl.de/).

```
iqtree -s $Indir/mafft.fasta -nt auto \
-pre $Outdir/out.fasta -mtree -m MFP -bb 1000 -bnni  -redo
```

| Arguments | Contents                                               |
| :-------- | ------------------------------------------------------ |
| -s        | File                                                   |
| -nt       | Number of threads                                      |
| -pre      | prefix of the resulting files                          |
| -mtree    | Input tree to infer site frequency model               |
| -m MFP    | Extended model selection followed by tree inference    |
| -bb       | Specify the fast method BS                             |
| -bnni     | Optimize UFBoot trees by NNI on bootstrap alignment    |
| -redo     | Ignore checkpoint and overwrite outputs (default: OFF) |

#### 10. iTOL

**Interactive Tree Of Life** is an online tool for the display and manipulation of phylogenetic trees. It provides most of the features available in other tree viewers, and offers a novel circular tree layout, which makes it easy to visualize mid-sized trees (up to several thousand leaves). Trees can be exported to several graphical formats, both bitmap and vector based.

![SARS-cov-2.png](https://github.com/Simple53/Markdown_Photos/blob/main/SARS-cov-2.png?raw=true) 

## **Questions**

#### **Pair-ended reads & Single-ended reads**

**Paired-ended Reads：**

- **Introduction**

    **Paired-end sequencing allows users to sequence both ends of a fragment and generate high-quality, alignable sequence data. Paired-end sequencing facilitates detection of genomic rearrangements and repetitive sequence elements, as well as gene fusions and novel transcripts.**
    
    In addition to producing twice the number of reads for the same time and effort in library preparation, sequences aligned as read pairs enable more accurate read alignment and the ability to detect insertion-deletion (indel) variants, which is not possible with single-read data. 

**![Pair-End reads](https://www.illumina.com/content/dam/illumina-marketing/images/science/v2/web-graphic/paired-end-vs-single-read-seq-web-graphic.jpg)**

- **Features**
  - **Simple Paired-End Libraries: Simple workflow allows generation of unique ranges of insert sizes**
   - **Efficient Sample Use: Requires the same amount of DNA as single-read genomic DNA or cDNA sequencing**
   - **Broad Range of Applications: Does not require methylation of DNA or restriction digestion; can be used for bisulfite sequencing**
   - **Simple Data Analysis: Enables high-quality sequence assemblies with short-insert libraries. A simple modification to the standard single-read library preparation process facilitates reading both the forward and reverse template strands of each cluster during one paired-end read. Both reads contain long-range positional information, allowing for highly precise alignment of reads.**

**Single-ended Reads:**

Single-read sequencing involves sequencing DNA from only one end, and is the simplest way to utilize Illumina sequencing. This solution delivers large volumes of high-quality data, rapidly and economically. Single-read sequencing can be a good choice for certain methods such as small RNA-Seq or chromatin immunoprecipitation sequencing (ChIP-Seq).   

**![See the source image](https://www.yourgenome.org/sites/default/files/illustrations/diagram/bioinformatics_single-end_pair-end_reads_yourgenome.png)**

#### **Sequencing Reads File** 

What is a FASTQ file？ 

**A FASTQ file usually contain millions of sequences. When these files are compressed with GZIP their sizes are reduced in more than 10 times (ZIP format is less efficient). Thus, the suffix will be ‘.fastq.gz’ or ‘.fq.gz’.**

**What does a FASTQ file look like?**

For each cluster that passes filter, a single sequence is written to the corresponding sample’s R1 FASTQ file. For a paired-end run, a single sequence is also written to the sample’s R2 FASTQ file. 

**![sequence](https://supportassets.illumina.com/content/dam/illumina-support/images/bulletins/11/fastq_files_explained_image.png)**

**Each entry in a FASTQ files consists of 4 lines:**

> 1. **Basic information. About the sequencing run and the cluster.**
>
>    ```
>    @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
>    ```
>
>    | **Element**            | **Requirements**                                     | **Description**                                              |
>    | ---------------------- | ---------------------------------------------------- | ------------------------------------------------------------ |
>    | **`@`**                | **@**                                                | **Each sequence  starts with @**                             |
>    | **`<instrument>`**     | **Characters allowed: a–z, A–Z, 0–9 and underscore** | **Instrument ID**                                            |
>    | **`<run number>`**     | **Numerical**                                        | **Run number on instrument**                                 |
>    | **`<flowcell ID>`**    | **Characters allowed: a–z, A–Z, 0–9**                |                                                              |
>    | **`<lane>`**           | **Numerical**                                        | **Lane number**                                              |
>    | **`<tile>`**           | **Numerical**                                        | **Tile number**                                              |
>    | **`<x_pos>`**          | **Numerical**                                        | **Run number on instrument**                                 |
>    | **`<y_pos>`**          | **Numerical**                                        | **X coordinate of cluster**                                  |
>    | **`<read>`**           | **Numerical**                                        | **Read number. 1 can be single read or Read 2 of paired-end** |
>    | **`<is filtered>`**    | **Y or N**                                           | **Y if the read is filtered (did not pass), N otherwise**    |
>    | **`<control number>`** | **Numerical**                                        | **0 when none of the control bits are on, otherwise it is an even number. On HiSeq X systems, control specification is not performed and this number is always 0.** |
>    | **`<sample number>`**  | **Numerical**                                        | **Sample number from sample sheet**                          |
>
> 2. **Sequence .(the base calls; A, C, T, G and N).**
>
> 3. **separator. A plus (+) sign.**
>
> 4. **Quality scores. These are Phred +33 encoded, using ASCII characters to represent the numerical quality scores.**

#### **Bash Scripts**

**Concept**

**A Bash script is a plain text file which contains a series of commands. These commands are a mixture of commands we would normally type ouselves on the command line. Anything you can run normally on the command line can be put into a script and it will do exactly the same thing. Similarly, anything you can put into a script can also be run normally on the command line and it will do exactly the same thing.**

Starting Off 

```
#!/bin/bash
A simple script
```

**The “#!” combo is called a shebang by most Unix geeks. This is used by the shell to decide which interpreter to run the rest of the script, and ignored by the shell that actually runs the script. Scripts can be written for all kinds of interpreters — bash, tsch, zsh, or other shells, or for Perl, Python, and so on.**

#### **Slurm job submission system**

**[Slurm](http://en.wikipedia.org/wiki/Simple_Linux_Utility_for_Resource_Management) is a resource manager and job scheduler designed to scheduled and allocated resources (CPU time, memory, etc.)**

**Gathering information**

```
 sinfo 
 # command gives an overview of the resources offered by the cluster
 # lists the partitions that are available
 # partition is a set of compute nodes
 # our default computing node is sugon 
 squeue
 # command shows to which jobs those resources are currently allocated.
```

**Creating a job**

  write a job script before submitting to the node. A job consists in two parts: resource requests and job steps:

- Resource requests: a number of CPUs, computing expected duration, amounts of RAM or disk space, etc. 

  > **`#!/bin/bash` for Linux to explain the script with bash language**
  > **`#SBATCH` directives just beneath the top command**
  > **`#SBATCH --job-name=test` job name**
  > **`#SBATCH --output=res.txt` log file**
  > **`#SBATCH --ntasks=1` task numbers**
  > **`#SBATCH --ntasks-per-node` task of each node**
  > **`#SBATCH --cpus-per-task` cpu of each task**
  > **`#SBATCH --mem-per-cpu=100` cpu memory**
  > **`#SBATCH --time=10:00` computing time**

- **Job step: tasks that must be done, software which must be run.**

#### BAM & SAM files

BAM and SAM formats are designed to contain the same information. The SAM format is more human readable, and easier to process by conventional text based processing programs, such as awk, sed, python, cut and so on. The BAM format provides binary versions of most of the same data, and is designed to compress reasonably well. 
