RefineM
=======

**[This project is in active development and not currently recommended for public use.]**


`version status <https://pypi.python.org/pypi/refinem>`_

`downloads <https://pypi.python.org/pypi/refinem>`_

RefineM is a set of tools for improving population genomes. It provides methods designed to improve the completeness of a genome along with methods for identifying and removing contamination. RefineM comprises only part of a full genome QC pipeline and should be used in conjunction with existing QC tools such as `CheckM <https://github.com/Ecogenomics/CheckM/wiki>`_. The functionality currently planned is:

*Improve completeness:*

* identify contigs with similarity to specific reference genome(s)

* identify contigs with compatible GC, coverage, and tetranucleotide signatures

* indetify partial population genomes which should be merged together (requires `CheckM <https://github.com/Ecogenomics/CheckM/wiki>`_)


*Reducing contamination:*

* taxonomically classify contigs within a genome in order to identify outliers

* identify contigs with divergent GC content, coverage, or tetranucleotide signatures

* identify contigs with a coding density suggestive of a Eukaryotic origin



Install
-------

The simplest way to install this package is through pip:
> sudo pip install refinem

This package requires numpy to be installed and makes use of the follow bioinformatic packages:

* `prodigal <http://prodigal.ornl.gov/>`_: Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. *Bioinformatics* 28: 2223-2230.

* `blast+ <http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_: Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. 2009. BLAST+: architecture and applications. *BMC Bioinformatics* 10:421: doi: 10.1186/1471-2105-10-421.

* `diamond <http://ab.inf.uni-tuebingen.de/software/diamond/>`_ Buchfink B, Xie C, Huson DH. 2015. Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12: 59–60 doi:10.1038/nmeth.3176.

* `krona <http://sourceforge.net/p/krona/home/krona/>`_ Ondov BD, Bergman NH, and Phillippy AM. 2011. Interactive metagenomic visualization in a Web browser. *BMC Bioinformatics* 12: 385.


Cite
----

If you find this package useful, please cite this git repository (https://github.com/dparks1134/refinem)


Copyright
---------

Copyright © 2015 Donovan Parks. See LICENSE for further details.
