## QIIME2 ITS protocol
# January 2023, Lady Grant
This is the workflow used for the Agbiome-Kerbel ITS dataset

Command cheat sheet: https://www.makeuseof.com/tag/mac-terminal-commands-cheat-sheet/

For Installation instructions: https://docs.qiime2.org

# Login instructions:
ssh zenith egrant@zenith.aggie.colostate.edu

password:O##M######K####

When you first log in you are in your directory, you will need to navigate to your porject directory.

The project directory is: Agribiome </home/projects/Agribiome>

All the original data coming fresh off the sequencer will be put here...


```python
/ORG-Data-wrighton/in_house_amplicon_libraries/...
```

# Activate QIIME environment on server prior to running anything
Once activated qiime2-2021.2 will show up before your name.

Ex: (qiime2-2021.2) [egrant@zenith Kerbel]$

If you don't activate the qiime environment the following error will be shown:

The error "bash: qiime: command not found"


```python
source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2
```

# 1. Import sequencing files into QIIME2
This step makes a directory named <reads> that will contain your barcodes and the forward and reverse primers.


```python
# mkdir <dir> = make directory < name of directory >
mkdir <reads> 
```

Copy the reads into your newly created reads folder and renames them to: <forward.fastq.gz>, <reverse.fastq.gz>, <barcodes.fastq.gz> .

Symbolic links from your reads into this directory will not use as much disk space as copying the reads. Both methods are shown below.

R1, R2 and I1 in the file name and correspond to read 1 (forward),  read 2 (reverse) and barcodes respectively.


```python
# symbolic link method...
ln -s [path to barcodes file] ./barcodes.fastq.gz
ln -s [path to forward reads file]  forward.fastq.gz
ln -s [path to reverse reads file]  reverse.fastq.gz
```


```python
# copy method
# cp <file> = copy file
# ./ <directory> = to here
# using the full file path is best here but not needed
# ex: < cp/file path to barcodes and primers ./reads/new name >
cp /ORG-Data-wrighton/in_house_amplicon_libraries/ITS_Kerbel_2021/Undetermined_S0_L001_R1_001.fastq.gz ./reads/forward.fastq.gz
cp /ORG-Data-wrighton/in_house_amplicon_libraries/ITS_Kerbel_2021/Undetermined_S0_L001_R2_001.fastq.gz ./reads/reverse.fastq.gz
cp /ORG-Data-wrighton/in_house_amplicon_libraries/ITS_Kerbel_2021/Undetermined_S0_L001_I1_001.fastq.gz ./reads/barcodes.fastq.gz
```

# 2. Demultiplex your sequences 
⏱ This step takes about 3 hours for this data set. Use screen and <&> so you can detach from screen later. 

link to screen commands: https://gist.github.com/jctosta/af918e1618682638aa82


Start a new session with new session name 


```python
# screen -S <session_name>
screen -S <session_name>
# activeate enviornement **new screen needs new environment** 
source /opt/Miniconda2/miniconda2/bin/activate qiime2-2021.2
```

Import Sequence files

If using EMP sequences which have one file for forward reads, one file for reverse reads and one file with associated barcodes and is still mulitplexed

For more info or if using other data format: https://docs.qiime2.org/2017.12/tutorials/importing/


```python
# using the & at the end of your commands will return your curser back!
# this command should take about 5 minutes
qiime tools import 
--type EMPPairedEndSequences 
--input-path reads/ 
--output-path emp-paired-end-sequences.qza &
```

Outout from command: emp-paired-end-sequences.qza

Before the next step you will need to upload your metadata file as a .txt file to the server. It must have the barcodes in a column labled "barcode-sequence" or at least match what is indicated in barcodes column flag.

** This step is different from the 16S protocol!**

Demultiplex the sequences


```python
qiime demux emp-paired 
--i-seqs emp-paired-end-sequences.qza # .qza file from above
--m-barcodes-file ITS_Metadata_qiime.txt #a .txt file that contains barcodes for the samples we want in the feature table.
--m-barcodes-column barcode-sequence # barcodes are in column labeled 'barcode-sequence'
--o-per-sample-sequences demux-ITS.qza  #output file containing reads separated by sample
--o-error-correction-details demux-details.qza 
--p-rev-comp-mapping-barcodes & #specifies barcodes need to be reverse complemented but varies based on sequencing center! If sequences per sample are low this can be adjusted.
```


```python
# Detatching screen, this will not kill the session.
# You can close your computer or work on other projects on the server.
screen -d <session_name>
```

Outputs from step: 

demux-ITS.qza

demux-details.qza

Visualize number of reads per sample and sequence quality.

After obtaining the .qzv file, you can visualize it in the QIIME2 view tools: https://view.qiime2.org/


```python
qiime demux summarize 
--i-data demux-ITS.qza # from previous step. Make sure these match.
--o-visualization demux-ITS.qzv
```

Outputs from step: demux-ITS.qzv  (interactive quality map)

# 3. Denoise the sequences (DADA2)
Denoises paired-end sequences, dereplicates them, and filters chimeras.

⏱ This step takes about 1.5 hours


```python
qiime dada2 denoise-single
--i-demultiplexed-seqs demux-ITS.qza #from previous step
--p-trunc-len 0
--p-max-ee 2
--p-trunc-q 2
--p-n-threads 20
--o-table table-ITS.qza
--o-representative-sequences rep-seqs-ITS.qza
--o-denoising-stats denoising-stats-ITS.qza
```

Outputs from step:

FeatureTable[Frequency]: table-ITS.qza

FeatureData[Sequence]: rep-seqs-ITS.qza 

SampleData[DADA2Stats]: denoising-stats-ITS.qza

Make and visualize feature table


```python
qiime feature-table summarize
  --i-table table-ITS.qza
  --o-visualization table-ITS.qzv
  --m-sample-metadata-file ITS_Metadata_qiime.txt #same metadata file uploaded eariler
```

Output from step: table-ITS.qzv

Make and visualize representative sequences


```python
qiime feature-table tabulate-seqs
--i-data rep-seqs-ITS.qza
--o-visualization rep-seqs-ITS.qzv
```

# 4. Training UNITE classifer
This step also takes a while. Use screen + <&>

Be sure to check that you are using the latest database. 

Downloads can be found here: https://unite.ut.ee/repository.php

Import UNITE data


```python
# run this in the <developer> folder downloaded from website
qiime tools import --type FeatureData[Sequence] 
--input-path sh_refs_qiime_ver9_99_s_29.11.2022_dev.fasta 
--output-path unite-seqs-ver9_99-11.2022.qza
```

Output from step: unite-seqs-ver9_99-11.2022.qza

Import taxonomy data


```python
qiime tools import --type FeatureData[Taxonomy] 
--input-path sh_taxonomy_qiime_ver9_97_s_29.11.2022_dev.txt 
--output-path unite-ver9-99-tax-11.2022.qza 
--input-format HeaderlessTSVTaxonomyFormat
```

Fit th classifer


```python
qiime feature-classifier fit-classifier-naive-bayes 
--i-reference-reads unite-seqs-ver9_99-11.2022.qza 
--i-reference-taxonomy unite-ver9-99-tax-11.2022.qza 
--o-classifier unite-ver9-99-classifier-11.2022.qza
```

Classify your reads


```python
qiime feature-classifier classify-sklearn 
--i-classifier unite-ver9-99-classifier-11.2022.qza
--i-reads rep-seqs-ITS.qza 
--o-classification taxonomy-ITS.qza
```

Visualize your taxonomy table


```python
qiime metadata tabulate
--m-input-file taxonomy-ITS.qza
--o-visualization taxonomy-ITS.qzv
```

Output from step: taxonomy-ITS.qzv

# Create ASV table
The next steps create a tab-delimited feature table that contains the read counts of each feature.


```python
qiime tools export 
--input-path taxonomy-ITS.qza 
--output-path feature_taxonomy_out
```


```python
qiime tools export 
--input-path table-ITS.qza 
--output-path feature_table_out
```


```python
cd feature_table_out
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv --header-key taxonomy
```

There are many different ways to combine taxonomy with ASV table but this will need to be done before any analysis in R.
