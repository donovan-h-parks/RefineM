###############################################################################
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import logging
import operator
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.common import (make_sure_path_exists,
                            alphanumeric_sort,
                            remove_extension,
                            concatenate_files)
from biolib.blast_parser import BlastParser
from biolib.external.diamond import Diamond
from biolib.taxonomy import Taxonomy
from biolib.plots.krona import Krona

from numpy import mean


"""
To Do:
 1. Need to get GC and coverage (and possibily tetras) into the mix
 2. Both the GeneProfile and Profile classes are very, very similar to classes
    in taxonomic_profile.py. These classes should be generalized and combined.
 *3. Need to enforce a taxonomically consistent classification for scaffolds in the
     taxa_profile and gene_profile methods. Without this, the Krona plot can give weird
     results (k__Bactera; p__Euryarchaeota). How to result this is an interesting question
     as the majority of genes could indicate Bacteria, but be spread over many different phyla
     so that Euryarchaeota is the majority classification at the phylum level.
 4. Report results should give profiles as 1) % scaffolds, 2) % genes, 3) % coding bases
 5. Make sure there is a table indicating the classification of each scaffold for use in identifying
    contamination
"""


class GeneProfile(object):
    """Create taxonomic profiles of genes across scaffolds within a genome.

    Genes are classified through homolog search against a
    database of reference genomes. Currently, homology search with
    diamond is supported. Each scaffold is assigned using a majority
    vote at each taxonomic rank. The taxonomic profile of a genome
    is given as the total number of scaffolds, total number of genes,
    and total number of coding bases assigned to each taxa.

    Note: this class deliberately ignores short scaffolds that do not
          contain at least one gene.
    """

    def __init__(self, cpus, output_dir):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        output_dir : str
            Directory to store results.
        """
        self.logger = logging.getLogger()

        self.cpus = cpus
        self.output_dir = output_dir

        self.rank_prefixes = Taxonomy().rank_prefixes
        self.rank_labels = Taxonomy().rank_labels

        # profile for each genome
        self.profiles = {}

    def taxonomic_profiles(self, table, taxonomy):
        """Create taxonomic profiles.

        Parameters
        ----------
        table : str
            Table containing hits to genomic fragments.
        taxonomy : dict[ref_genome_id] -> [domain, phylum, ..., species]
            Taxonomic assignment of each reference genome.
        """

        blast_parser = BlastParser()

        hit_summary = defaultdict(lambda: defaultdict(int))
        for hit in blast_parser.read_hit(table):
            seq_id, genome_id = hit.query_id.split('~')
            scaffold_id = seq_id[0:seq_id.rfind('_')]
            ref_genome_id = hit.subject_id[hit.subject_id.rfind('~') + 1:]

            self.profiles[genome_id].add_hit(scaffold_id,
                                             taxonomy[ref_genome_id],
                                             hit.evalue,
                                             hit.perc_identity,
                                             hit.aln_length)

            hit_summary[genome_id][ref_genome_id] += 1

    def write_genome_summary(self, output_file):
        """Summarize classification of each genome.

        Parameters
        ----------
        output_file : str
            Output file.
        """

        fout = open(output_file, 'w')
        fout.write('Genome id\tLength (bp)\t# sequences')
        for rank in self.rank_labels:
            fout.write('\t' + rank + ': taxa')
            fout.write('\t' + rank + ': percent of bps')
            fout.write('\t' + rank + ': percent of sequences')
            fout.write('\t' + rank + ': avg. evalue')
            fout.write('\t' + rank + ': avg. perc identity')
            fout.write('\t' + rank + ': avg. align length (AA)')
        fout.write('\n')

        sorted_genome_ids = alphanumeric_sort(self.profiles.keys())
        for genome_id in sorted_genome_ids:
            self.profiles[genome_id].write_genome_summary(fout)

        fout.close()

    def run(self, gene_files, db_file, taxonomy_file, evalue, per_identity):
        """Create taxonomic profiles for a set of genomes.

        Parameters
        ----------
        gene_files : list of str
            Fasta files of called genes to process.
        db_file : str
            Database of reference genes.
        taxonomy_file : str
            File containing GreenGenes taxonomy strings for reference genomes.
        evalue : float
            E-value threshold used by blast.
        per_identity: float
            Percent identity threshold used by blast.
        """

        # read taxonomy file
        self.logger.info('')
        self.logger.info('  Reading taxonomic assignment of reference genomes.')
        taxonomy = Taxonomy().read(taxonomy_file)

        # concatenate gene files
        diamond_output_dir = os.path.join(self.output_dir, 'diamond')
        make_sure_path_exists(diamond_output_dir)

        gene_file = os.path.join(diamond_output_dir, 'genes.faa')
        concatenate_files(gene_files, gene_file)

        # record length and number of genes in each scaffold
        for aa_file in gene_files:
            genome_id = remove_extension(aa_file, '.genes.faa')
            self.profiles[genome_id] = Profile(genome_id)

            for seq_id, seq in seq_io.read_seq(aa_file):
                seq_id = seq_id[0:seq_id.rfind('~')]
                scaffold_id = seq_id[0:seq_id.rfind('_')]
                self.profiles[genome_id].genes_in_scaffold[scaffold_id] += 1
                self.profiles[genome_id].coding_bases[scaffold_id] += len(seq)

        # run diamond and create taxonomic profile for each genome
        self.logger.info('')
        self.logger.info('  Running diamond blastp with %d processes (be patient!)' % self.cpus)

        diamond = Diamond(self.cpus)
        diamond_daa_out = os.path.join(diamond_output_dir, 'diamond_hits')
        diamond.blastp(gene_file, db_file, evalue, per_identity, 1, diamond_daa_out)

        diamond_table_out = os.path.join(diamond_output_dir, 'diamond_hits.tsv')
        diamond.view(diamond_daa_out + '.daa', diamond_table_out)

        # create taxonomic profile for each genome
        self.logger.info('')
        self.logger.info('  Creating taxonomic profile for each genome.')
        self.taxonomic_profiles(diamond_table_out, taxonomy)

        # write out taxonomic profile
        self.logger.info('')
        self.logger.info('  Writing taxonomic profile for each genome.')
        report_dir = os.path.join(self.output_dir, 'bin_reports')
        make_sure_path_exists(report_dir)
        for genome_id, profile in self.profiles.iteritems():
            seq_summary_out = os.path.join(report_dir, genome_id + '.sequences.tsv')
            profile.write_seq_summary(seq_summary_out)

            genome_profile_out = os.path.join(report_dir, genome_id + '.profile.tsv')
            profile.write_genome_profile(genome_profile_out)

        # create summary report for all genomes
        genome_summary_out = os.path.join(self.output_dir, 'genome_summary.tsv')
        self.write_genome_summary(genome_summary_out)

        # create Krona plot
        krona_profiles = defaultdict(lambda: defaultdict(int))
        for genome_id, profile in self.profiles.iteritems():
            seq_assignments = profile.classify_seqs()

            for seq_id, classification in seq_assignments.iteritems():
                taxa = []
                for r in xrange(0, len(profile.rank_labels)):
                    taxa.append(classification[r][0])

                krona_profiles[genome_id][';'.join(taxa)] += profile.genes_in_scaffold[seq_id]

        krona = Krona()
        krona_output_file = os.path.join(self.output_dir, 'taxonomic_profiles.krona.html')
        krona.create(krona_profiles, krona_output_file)


class Profile(object):
    """Profile of hits to reference genomes."""

    def __init__(self, genome_id):
        """Initialization."""

        self.percent_to_classify = 0.2

        self.rank_prefixes = Taxonomy().rank_prefixes
        self.rank_labels = Taxonomy().rank_labels

        self.genome_id = genome_id
        self.unclassified = 'unclassified'

        self.TaxaInfo = namedtuple('TaxaInfo', """evalue
                                                perc_identity
                                                aln_length
                                                num_seqs
                                                num_basepairs""")

        # track hits at each rank: dict[contig_id][rank][taxa] -> [HitInfo, ...]
        self.HitInfo = namedtuple('HitInfo', """evalue
                                                perc_identity
                                                aln_length""")
        self.hits = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        # number of coding bases in scaffold
        self.coding_bases = defaultdict(int)

        # number of fragments from each sequence
        self.genes_in_scaffold = defaultdict(int)

    def add_hit(self, seq_id, taxonomy, evalue, perc_identity, aln_length):
        """Add hit to profile.

        Parameters
        ----------
        seq_id : str
            Unique identifier of sequence from which hit arose.
        taxonomy : list
            List indication taxa at each taxonomic rank.
        evalue : float
            E-value of hit.
        per_identity: float
            Percent identity of hit.
        aln_length: int
            Alignment length of hit.
        """

        d = self.hits[seq_id]
        for i, taxa in enumerate(taxonomy):
            d[i][taxa].append(self.HitInfo(evalue, perc_identity, aln_length))

    def classify_seqs(self):
        """Classify sequences.

        Sequences are classified using a majority vote
        over all fragments originating from the sequence
        with a valid hit. If less than 20% of fragments have
        a valid hit, the sequence is considered unclassified.
        Classification is performed from the highest (domain)
        to lowest (species) rank. If a rank is taxonomically
        inconsistent with a higher ranks classification, this
        rank and all lower ranks are set to unclassified.

        Returns
        -------
        dict : d[contig_id][rank] -> [taxa, HitInfo]
            Classification of each sequence along with summary statistics
            of hits to the specified taxa.
        """

        # classify each sequence using a majority vote
        seq_assignments = defaultdict(lambda: defaultdict(list))
        for seq_id, rank_hits in self.hits.iteritems():
            for rank in xrange(0, len(self.rank_prefixes)):
                taxa = max(rank_hits[rank], key=lambda x: len(rank_hits[rank][x]))
                count = len(rank_hits[rank][taxa])

                if count >= self.percent_to_classify * self.genes_in_scaffold[seq_id]:
                    seq_assignments[seq_id][rank] = [taxa, rank_hits[rank][taxa]]
                else:
                    seq_assignments[seq_id][rank] = [self.unclassified, None]

        # identify sequences with no hits
        for seq_id in self.genes_in_scaffold:
            if seq_id not in seq_assignments:
                for rank in xrange(0, len(self.rank_prefixes)):
                    seq_assignments[seq_id][rank] = [self.unclassified, None]

        return seq_assignments

    def profile(self):
        """Relative abundance profile at each taxonomic rank.

        Relative abundance is derived from the number
        of base pairs assigned to a given taxa.

        Parameters
        ----------
        rank : int
            Desired rank.

        Returns
        -------
        dict : d[rank][taxa] -> percentage
            Relative abundance of taxa at a given rank.
        dict : d[rank][taxa] -> TaxaInfo
           Statistics for each taxa.
        """

        seq_assignments = self.classify_seqs()

        total_coding_bps = sum(self.coding_bases.values())

        profile = defaultdict(lambda: defaultdict(float))
        stats = defaultdict(dict)

        for r in xrange(0, len(self.rank_labels)):
            num_seqs = defaultdict(int)
            num_basepairs = defaultdict(int)
            hit_stats = defaultdict(list)
            for seq_id, data in seq_assignments.iteritems():
                taxa, hit_info = data[r]
                profile[r][taxa] += float(self.coding_bases[seq_id]) / total_coding_bps
                num_seqs[taxa] += 1
                num_basepairs[taxa] += self.coding_bases[seq_id]

                if taxa != self.unclassified:
                    hit_stats[taxa].extend(hit_info)

            # calculate averages of hit statistics
            for taxa, hit_info in hit_stats.iteritems():
                avg_evalue = mean([x.evalue for x in hit_info])
                avg_perc_identity = mean([x.perc_identity for x in hit_info])
                avg_aln_length = mean([x.aln_length for x in hit_info])

                stats[r][taxa] = self.TaxaInfo(avg_evalue,
                                            avg_perc_identity,
                                            avg_aln_length,
                                            num_seqs[taxa],
                                            num_basepairs[taxa])

            stats[r][self.unclassified] = self.TaxaInfo(None,
                                                     None,
                                                     None,
                                                     num_seqs[self.unclassified],
                                                     num_basepairs[self.unclassified])

        return profile, stats

    def write_genome_summary(self, fout):
        """Write most abundant taxa at each rank.

        Parameters
        ----------
        fout : output stream
            Output stream.
        """

        profile, stats = self.profile()

        fout.write('%s\t%d\t%d' % (self.genome_id, sum(self.coding_bases.values()), len(self.coding_bases)))
        for r in xrange(0, len(self.rank_labels)):
            taxa, percent = max(profile[r].iteritems(), key=operator.itemgetter(1))

            total_seq = sum([stats[r][t].num_seqs for t in stats[r]])

            if taxa != self.unclassified:
                fout.write('\t%s\t%.2f\t%.2f\t%.1g\t%.1f\t%.1f' % (taxa,
                                                              percent * 100,
                                                              stats[r][taxa].num_seqs * 100.0 / total_seq,
                                                              stats[r][taxa].evalue,
                                                              stats[r][taxa].perc_identity,
                                                              stats[r][taxa].aln_length))
            else:
                fout.write('\t%s\t%.2f\t%.2f\t%s\t%s\t%s' % (taxa,
                                                          percent * 100,
                                                          stats[r][taxa].num_seqs * 100.0 / total_seq,
                                                          'na',
                                                          'na',
                                                          'na'))

        fout.write('\n')

    def write_genome_profile(self, output_file):
        """Write complete profile for genome.

        Parameters
        ----------
        output_file : str
            Output file.
        """

        fout = open(output_file, 'w')
        for rank in self.rank_labels:
            if rank != self.rank_labels[0]:
                fout.write('\t')

            fout.write(rank + ': taxa')
            fout.write('\t' + rank + ': percent of coding bps')
            fout.write('\t' + rank + ': percent of sequences')
            fout.write('\t' + rank + ': avg. evalue')
            fout.write('\t' + rank + ': avg. perc identity')
            fout.write('\t' + rank + ': avg. align length (aa)')
        fout.write('\n')

        # sort profiles in descending order of abundance
        profile, stats = self.profile()

        sorted_profiles = {}
        max_taxa = 0
        for r in xrange(0, len(self.rank_labels)):
            sorted_profile = sorted(profile[r].items(), key=operator.itemgetter(1))
            sorted_profile.reverse()

            sorted_profiles[r] = sorted_profile

            if len(sorted_profiles) > max_taxa:
                max_taxa = len(sorted_profiles)

        # write out table
        for i in xrange(0, max_taxa):
            for r in xrange(0, len(self.rank_labels)):
                if r != 0:
                    fout.write('\t')

                if len(sorted_profiles[r]) > i:
                    total_seq = sum([stats[r][t].num_seqs for t in stats[r]])

                    taxa, percent = sorted_profiles[r][i]

                    if taxa != self.unclassified:
                        fout.write('%s\t%.2f\t%.2f\t%.1g\t%.1f\t%.1f' % (taxa,
                                                                      percent * 100,
                                                                      stats[r][taxa].num_seqs * 100.0 / total_seq,
                                                                      stats[r][taxa].evalue,
                                                                      stats[r][taxa].perc_identity,
                                                                      stats[r][taxa].aln_length))
                    else:
                        fout.write('%s\t%.2f\t%.2f\t%s\t%s\t%s' % (taxa,
                                                                  percent * 100,
                                                                  stats[r][taxa].num_seqs * 100.0 / total_seq,
                                                                  'na',
                                                                  'na',
                                                                  'na'))
                else:
                    fout.write('\t\t\t\t\t')

            fout.write('\n')

    def write_seq_summary(self, output_file):
        """Summarize classification of each sequence.

        Parameters
        ----------
        output_file : str
            Output file.
        """

        seq_assignments = self.classify_seqs()

        fout = open(output_file, 'w')
        fout.write('Sequence id\tCoding bases (bp)\t# genes')
        for rank in self.rank_labels:
            fout.write('\t' + rank + ': taxa')
            fout.write('\t' + rank + ': hits (%)')
            fout.write('\t' + rank + ': avg. evalue')
            fout.write('\t' + rank + ': avg. perc identity')
            fout.write('\t' + rank + ': avg. align length (AA)')
        fout.write('\n')

        for seq_id in seq_assignments:
            fout.write('%s\t%d\t%d' % (seq_id,
                                       self.coding_bases[seq_id],
                                       self.genes_in_scaffold[seq_id]))

            for r in xrange(0, len(self.rank_labels)):
                taxa, hit_info = seq_assignments[seq_id][r]

                if taxa != self.unclassified:
                    avg_evalue = mean([x.evalue for x in hit_info])
                    avg_perc_identity = mean([x.perc_identity for x in hit_info])
                    avg_aln_length = mean([x.aln_length for x in hit_info])

                    hit_str = '%d (%.1f%%)' % (len(hit_info),
                                               len(hit_info) * 100.0 / self.genes_in_scaffold[seq_id])
                    fout.write('\t%s\t%s\t%.1g\t%.1f\t%.1f' % (taxa,
                                                               hit_str,
                                                               avg_evalue,
                                                               avg_perc_identity,
                                                               avg_aln_length))
                else:
                    fout.write('\t%s\t%s\t%s\t%s\t%s' % (taxa,
                                                           'na',
                                                           'na',
                                                           'na',
                                                           'na'))
            fout.write('\n')
