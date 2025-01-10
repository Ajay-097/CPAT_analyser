# CPAT lncRNA analyser

CPAT lncRNA analyser is a mini pipeline used to analyse the results from the CPAT tool (Coding Potential Assessment Tool) with the context of identifying long non-coding RNAs. The main inputs required for this pipeline is a .tsv file containanig the probabilities for each of the ORFs (Open reading frames) identified using CPAT and a .log file that has the transcript IDs for which no ORFs have been identified.
The tool produces a .txt file with the transcript IDs of all the potential lncRNAs identified using CPAT.


The celegans_cpat.ORF_prob.tsv and celegans_CPAT_run_info.log here are generated for _Caenorhabditis elegans_ as deltailed below. Using the two files the CPAT lncRNA analyser pipeline gives the transcript IDs of the candidate lncRNAs which can be used to annotate lncRNAs. 


## How to obtain the input files?
Here we detail the steps on how to run CPAT and obtain the prob.tsv file and the non ORF transcript IDs. Before we run CPAT we perform some manual filtering steps to remove any transcript that are obviously not lncRNAs.
### 1. Filtering mono exonic transcripts, Getting a .txt file with the transcript Ids that have more than 2 exons
We use a tool such as [StringTie2](https://github.com/skovaka/stringtie2) to assemble transcripts from RNA seq data.
```
awk '$3=="exon"' stringtie2_transcripts.gtf | \
awk '{match($0, /transcript_id "([^"]+)"/, a); print a[1]}' | \
sort | uniq -c | awk '$1 >=2 {print $2}' > more_than_2_exons.txt
```
### 2. Retaining only transcripts>= 2 exons in the gtf using the .txt file from the previous step
```
# Use the -w flag for when you need only the exact matches, otherwise we might get partial matching
grep -F -w -f <(sed 's/^/transcript_id "/; s/$/";/' more_than_2_exons.txt) stringtie2_transcripts.gtf > transcripts_exon_filtered.gtf
```
### 3. Filtering the transcripts overlapping with annotated exons
```
# Convert GTF files to BED for overlap comparison
awk '$3 == "exon" {print $1, $4, $5, $9}' OFS="\t" transcripts_exon_filtered.gtf > transcripts_exons.bed
awk '$3 == "exon" {print $1, $4, $5, $9}' OFS="\t" reference_annotations.gtf > annotated_exons.bed

# Use bedtools to remove overlapping transcripts
bedtools intersect -v -a transcripts_exons.bed -b annotated_exons.bed | cut -f4 > transcripts_no_overlap.txt

# The above .txt will have duplicates and other characters like " and ; in the transcript_id
awk '{match($0, /"([^"]+)"/, a); print a[1]}' transcripts_no_overlap.txt| sort | uniq > temp.txt
rm transccripts_no_overlap.txt
mv temp.txt transcripts_no_overlap.txt

# Filter out the overlapping exons from the .gtf from step 2 using the transccripts_no_overlap.txt file 
grep -F -w -f <(sed 's/^/transcript_id "/; s/$/";/' transcripts_no_overlap.txt) transcripts_exon_filtered.gtf> transcripts_exon_filtered_no_overlap.gtf
```

### 4. Running CPAT
- Download the CPAT tool based on instructions from the [CPAT manual](https://cpat.readthedocs.io/en/latest/#run-cpat-on-local-computer).
- We need a hexamer table and a logistic regression model set up as inputs to run CPAT. They have prebuilt this for model species but for any other species, we need to run the scripts provided with the CPAT tool to generate these.
- Get the fasta sequences from the transcripts_exon_filtered_no_overlap.gtf using any tool such as [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread_ex). 

```
# Creating a hexamer table using a fasta inputs of coding and non-coding sequences. This is to train the model to differentiate certain hexamers as coding and non-coding

make_hexamer_tab -c coding_sequences.fa -n noncoding_sequences.fa > sratti_hexamer.tsv

# Creating the logit
# We need to supply a fasta input of mRNA transcripts and non-coding transcripts as input along with the hexamer table.
make_logitModel -x sratti_hexamer.tsv -c <mRNA.fasta> -n <non_coding_transcripts.fasta> -o <output_prefix>

# Running CPAT
cpat -x Hexamer.tsv  -d  logitModel.RData  --top-orf=100  --antisense -g transcripts_exon_filtered_no_overlap.fa -o output
```
### NOTE
- You must specify `-antisense`, otherwise, it will only search ORFs from the sense strand.
- You also specify `-top-orf` to a large number to report all the ORFs.
- The `-min-orf` is set to 75 by default, same as [NCBI ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/).

