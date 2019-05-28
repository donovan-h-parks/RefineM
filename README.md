# RefineM

**[This project is in active development. Documentation is fairly light. You are welcomed to use this software, but please expect it to change in non-trivial ways.]**

_All users are encouraged to update to v0.0.22. In previous versions, both a mean absolute error and correlation criteria were used to identify contigs with divergent coverage profiles. Starting with v0.0.21, only the mean absolute error criteria is used by default. The correlation criteria can be misleading with fewer than 6 data points (i.e., BAM files) so is not used by default. Use of the ssu_erroneous method requires v0.0.22._

[![version status](https://img.shields.io/pypi/v/refinem.svg)](https://pypi.python.org/pypi/refinem)

RefineM is a set of tools for improving population genomes. It provides methods designed to improve the completeness of a genome along with methods for identifying and removing contamination. RefineM comprises only part of a full genome QC pipeline and should be used in conjunction with existing QC tools such as [CheckM](https://github.com/Ecogenomics/CheckM/wiki). The functionality currently planned is:

*Reducing contamination:*
* taxonomically classify contigs within a genome in order to identify outliers
* identify contigs with divergent GC content, coverage, or tetranucleotide signatures
* identify contigs with a coding density suggestive of a Eukaryotic origin (in progress)

*Improve completeness (in progress):*
* identify contigs with similarity to specific reference genome(s)
* identify contigs with compatible GC, coverage, and tetranucleotide signatures
* identify partial MAGs which should be merged together (requires [CheckM](https://github.com/Ecogenomics/CheckM/wiki))

## Install

The simplest way to install this package is through pip:
> sudo pip install refinem

This package requires numpy to be installed and makes use of the follow bioinformatic packages:

* [prodigal](http://prodigal.ornl.gov/) >=2.6.3: Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. *Bioinformatics* 28: 2223-2230.
* [blast+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) >=2.6.0: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. *BMC Bioinformatics* 10:421: doi: 10.1186/1471-2105-10-421.
* [diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) >=0.9.9: Buchfink B, Xie C, Huson DH. 2015. Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12: 59–60 doi:10.1038/nmeth.3176.
* [krona](http://sourceforge.net/p/krona/home/krona/) >=2.7: Ondov BD, Bergman NH, and Phillippy AM. 2011. Interactive metagenomic visualization in a Web browser. *BMC Bioinformatics* 12: 385.

## Identifying potential contamination

RefineM can identify potential contamination based on the genomic properties (GC, tetranucleotide signatures, coverage) of scaffolds and based on their taxonomic assignment against a reference database.

### Removing contamination based on genomic properties

To identify scaffolds with genomic properties that are divergent from the expect values for a bin (i.e., metagenome-assembled genome or MAG), the tetranucleotide signature and coverage profiles for scaffolds must be calculated:
```
>refinem scaffold_stats -c 16 <scaffold_file> <bin_dir> <stats_output_dir> <bam_files>
```
where <scaffold_file> is a FASTA file containing the scaffolds/contigs binned to produce your bins, <bin_dir> is the directory containing your bins, <stats_output_dir> is the directory to store results, and <bam_files> is one or more indexed BAM files specifying the mapping of reads to the scaffolds/contigs. The number of CPUs to use can be specified with the -c flag.

Scaffolds with divergent genomic properties can then be identified using:
```
>refinem outliers <stats_output_dir>/scaffold_stats.tsv <outlier_output_dir>
```
where the scaffold_stats.tsv file is produced by the scaffold_stat command and the <outlier_output_dir> will contain a number of data files and plots for manually investigating the genomic properties of scaffolds within your bins. Starting with RefineM v0.0.21, contigs with divergent coverage profiles are identified exclusively using a mean absolute percent error criteria (--cov_perc). If you have more than 6 data point (i.e. BAM files) comprising your coverage profiles you may wish to consider using the coverage correlation criteria (--cov_corr) instead of or in addition to this absolute error criteria. The output file outliers.tsv indicates the scaffolds RefineM has identified as being potential contamination. If desired, you can modify the criteria used by RefineM to identify potential contamination (see refinem outliers -h).

Contaminating scaffolds can be removed from your bins as follows:
```
>refinem filter_bins <bin_dir> <outlier_output_dir>/outliers.tsv <filtered_output_dir>
```
where <bin_dir> is the directory containing your bins to be modified, outliers.tsv indicates the scaffolds to remove from each bin and is produced by the outliers command, and <filtered_output_dir> will contain your bins with the specified scaffolds removed. If your only want the output directory to contain bins that were modified, you can use the --modified_only flag.

### Removing contamination based on taxonomic assignments

To identify scaffolds with taxonomic assignments that are divergent from the taxonomic affliations of a bin, the genes in each scaffold/contig are classified against a reference database using DIAMOND. You can call genes on all your genomes using:
```
>refinem call_genes -c 40 <bin_dir> <gene_output_dir>
```
where <bin_dir> is the directory containing your bins and the directory <gene_output_dir> will contain called genes for your bins. 

The genes comprising each bin can then be classified against a reference database using:
```
>refinem taxon_profile -c 40 <gene_output_dir> <stats_output_dir>/scaffold_stats.tsv <reference_db> <reference_taxonomy> <taxon_profile_output_dir>
```
where <gene_output_dir> is the output of the call_genes command, <stats_output_dir>/scaffold_stats.tsv is the output from the scaffold_stats command as discussed [above](#removing-contamination-based-on-genomic-properties), and the <reference_db> and <reference_taxonomy> are used as reference database for assigning taxonomic classifications to individual genes based on a top hit criteria. Reference files and their format are discussed [below](#reference-database-and-taxonomy-files). 

Scaffolds with divergent taxonomic assignments can then be identified with:
```
>refinem taxon_filter -c 40 <taxon_profile_dir> taxon_filter.tsv
```
where <taxon_profile_dir> is the output directory of the taxon_profile command and scaffolds determined to be contamination are written to taxon_filter.tsv. If desired, you can modify the criteria used by RefineM to identify potential contamination (see refinem taxon_filter -h).

Contaminating scaffolds can be removed from your bins as follows:
```
>refinem filter_bins <bin_dir> taxon_filter.tsv <filtered_output_dir>
```
where <bin_dir> is the directory containing your bins to be modified, taxon_filter.tsv indicates the scaffolds to remove from each bin and is produced by the taxon_filter command, and <filtered_output_dir> will contain your bins with the specified scaffolds removed. If your only want the output directory to contain bins that were modified, you can use the --modified_only flag.

### Removing contigs with incongruent 16S rRNA genes (requires version >=0.0.22)

Scaffolds with 16S rRNA genes that appear incongruent with the taxonomic identity of a bin can be identified as follows:
```
>refinem ssu_erroneous <bin_dir> <taxon_profile_dir> <ssu_db> <reference_taxonomy> <ssu_output_dir>
```
where <bin_dir> is the directory containing your bins, <taxon_profile_dir> is the output directory of the taxon_profile command, and the <ssu_db> and <reference_taxonomy> are reference database for establishing the taxonomic identity of 16S rRNA genes. Reference files and their format are discussed [below](#reference-database-and-taxonomy-files). Output files will be placed in the <ssu_output_dir>. The file ssu_erroneous.tsv lists genomes and corresponding scaffolds that may have erroneous 16S rRNA genes. The genome classification indicates the percentage of genes which support each taxon assignment (e.g., d__Bacteria (99%)) and the 16S rRNA classification indicates the taxonomic assignment of the top hit in the reference 16S rRNA database. The E-value, alignment length, and percent identity of the top hit is provide in order to allow the quality of the top hits to be assessed. It is recommended that this file be inspected and your judgement be used to decide which contigs should be deemed erroneous. Contigs specified in the ssu_erroneous.tsv file can be removed from bins using the filter_bins method of RefineM.

## Reference database and taxonomy files

Reference protein and 16S rRNA databases and a corresponding taxonomy file can be obtained from https://data.ace.uq.edu.au/public/misc_downloads/refinem/. These databases are constructed from the dereplicated set of genomes used to define the Genome Taxonomy Database (GTDB: http://gtdb.ecogenomic.org/). The taxonomy file can be used with both the protein and 16S rRNA databases.

If you wish to make you own reference protein database, the protein sequences must be formatted into a DIAMOND database and have header information in the format <genome_id>~<contig_id>_<gene_num>, e.g.:
```
>GCF_001687105.1~contig000001_1
```

If you wish to make you own 16S rRNA database, the sequences must be formatted into a BLASTN database and have the genome identifier in the FASTA header, e.g.:
```
>GCF_001687105.1
``` 

Taxonomy information for reference genomes must be provided in a seperate taxonomy file. This file is a simple tab-separated values file with two columns indicating the genome ID and Greengenes-style taxonomy string, e.g.:
```
>GCF_001687105.1    d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhodobacterales;f__Rhodobacteraceae;g__Yangia;s__
```

## Cite

If you find this package useful, please cite:

Parks DH et al. 2017. [Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life](http://dx.doi.org/10.1038/s41564-017-0012-7). Nat Microbiol 2: 1533-1542.

Please also consider citing the 3rd party applications required by RefineM such as [Prodigal](http://prodigal.ornl.gov/) and [DIAMOND](http://ab.inf.uni-tuebingen.de/software/diamond/).

## Copyright

Copyright © 2015 Donovan Parks. See LICENSE for further details.
