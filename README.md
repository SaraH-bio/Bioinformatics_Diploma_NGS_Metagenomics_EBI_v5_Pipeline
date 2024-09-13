**Bioinformatics_Diploma_NGS_Metagenomics_EBI_v5_Pipeline**

# download silva dbs
mkdir silva_ssu silva_lsu
wget \
  ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_ssu-20200130.tar.gz \
  ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/silva_lsu-20200130.tar.gz 
tar --extract --gzip --directory=silva_ssu silva_ssu-20200130.tar.gz
tar --extract --gzip --directory=silva_lsu silva_lsu-20200130.tar.gz

# download Pfam ribosomal models
mkdir ribosomal
wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/rfam_models/ribosomal_models/ribo.cm \
  ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs/rfam_models/ribosomal_models/ribo.claninfo \
  -P ribosomal 

# Download SeqPrep
wget https://github.com/jstjohn/SeqPrep/archive/refs/heads/master.zip -O SeqPrep-master.zip
unzip SeqPrep-master.zip
cd SeqPrep-master
make

# Download sratoolkit
conda install daler::sratoolkit


**Pipeline Steps**
# Download the Sample Using Prefetch
prefetch SRR5903755

# Split Sample into FASTQ Files Using Fastq-dump
fastq-dump SRR5903755 --split-files

# Merge Reads Using SeqPrep
for R1 in *1.fastq; do
    R2="${R1%_1.fastq}_2.fastq"
    SeqPrep -f "$R1" -r "$R2" -1 "unmerged_${R1}" -2 "unmerged_${R2}" -s "${R1%_1.fastq}_fastq_merged.gz"
done

# FastQC before Trimming step (Quality check) 
for i in *merged.gz; do
    fastqc "$i"
done

# Trim the merged files using Trimmomatic
for i in *merged.gz; do
    trimmomatic SE -phred33 "$i" "${i%.gz}_trimmed.fq.gz" ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100
done

# FastQC after Trimming step (Quality check) 
for i in *trimmed.fq.gz; do
    fastqc "$i"
done

# Unzip trimmed files
gunzip ./*trimmed.fq.gz

 # Convert FASTQ to FASTA
for i in *trimmed.fq; do
    seqkit fq2fa "$i" > "${i%.fq}.fa"
done

# Create directory for results
mkdir -p cmresults

# Search for homologous sequences using CMsearch
for i in ribsome/*.cm; do
    base=$(basename "$i" .cm)
    for trimmed in *trimmed.fa; do
        trimmed_base=$(basename "$trimmed" .fa)
        cmsearch --tblout "./cmresults/${base}_${trimmed_base}.tbl" "$i" "$trimmed"
    done
done

 # Process tbl files for sample 1
for tbl_file in cmresults/*SRR5903755_fastq_merged_trimmed.tbl; do
    awk '/^[^#]/ {if ($8 > $9) {print $1 "\t" $9-1 "\t" $8} else {print $1 "\t" $8-1 "\t" $9}}' "$tbl_file" >>combined_SRR5903755.bed
done

# Mask FASTA files using BED files
bedtools maskfasta -fi SRR5903755_fastq_merged_trimmed.fa -bed combined_SRR5903755.bed -fo 3755_masked.fast
a

# Extract FASTA sequences from BED files
bedtools getfasta -fi SRR5903755_fastq_merged_trimmed.fa -bed combined_SRR5903755.bed -fo 3755_get.fasta

# Map sequences to Silva LSU and SSU databases
# Map sequences to LSU database
mapseq 3755_get.fasta silva_lsu-20200130/LSU.fasta slv_lsu_filtered2.txt > 3755lsu.txt

# Map sequences to SSU database
mapseq 3755_get.fasta silva_ssu-20200130/SSU.fasta slv_ssu_filtered2.txt > 3755ssu.txt

# Count OTUs and convert to Krona format
for i in *su.txt; do
    mapseq -otucounts "$i" > "${i%.txt}.otu"
done

for i in *.otu; do
    awk 'NR>1 {gsub(/;/, "\t", $3); print $4 "\t" $3}' "$i" > "${i%.otu}krona.txt"
done

for i in *krona.txt; do
    ktImportText "$i" -o "${i%.txt}.html"
done
