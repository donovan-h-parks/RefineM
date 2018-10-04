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
                            remove_extension)
from biolib.blast_parser import BlastParser
from biolib.external.execute import check_dependencies
from biolib.external.diamond import Diamond
from biolib.taxonomy import Taxonomy
from biolib.plots.krona import Krona

from refinem import version
from refinem.common import concatenate_gene_files
from refinem.scaffold_stats import ScaffoldStats


"""
To Do:
    Get some form of gene annotation into the mix?
"""

def mean(v):
    """Calculate mean of values.
    
    This is a weird workaround. Ideally, numpy.mean would be used. However,
    this is causing an issue with small numbers, e.g. mean([0.0, 5.3e-10]).
    This seems to work if entered directly into the Python console, but 
    not while running RefineM. I have been unable to determine why.
    """
    
    return sum(v) / len(v)


class TaxonProfile(object):
    """Create taxonomic profiles of genes across scaffolds within a genome.

    Genes are classified through homolog search against a
    database of reference genomes. Currently, homology search with
    diamond is supported. Each scaffold is classified as follows:

    1. starting at the rank of domain, the scaffold is assigned to
        the taxon with the most votes unless fewer than X% (user specific)
        of the genes are assigned to a single taxon. In this case, the
        scaffold is marker as unclassified.
    2. Lower ranks are assigned in the same manner, except that if the
        assigned taxon in not taxonomically consistent with previously
        assigned taxon than this rank and all lower ranks are set to
        unclassified.

    The taxonomic profile of a genome is given as:

    1. the total number of scaffolds assigned to a taxon, weighted
        by the number of genes in each scaffold
    2. the total number of genes assigned to a taxon without
        regard to the source scaffold of each gene
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

        self.logger = logging.getLogger('timestamp')
        
        check_dependencies(('diamond', 'ktImportText'))

        self.cpus = cpus
        self.output_dir = output_dir

        # profile for each genome
        self.profiles = {}
        
    def taxonomic_profiles(self, table, taxonomy):
        """Create taxonomic profiles.

        Parameters
        ----------
        table : str
            Table containing hits to genes.
        taxonomy : d[ref_genome_id] -> [domain, phylum, ..., species]
            Taxonomic assignment of each reference genome.
        """

        blast_parser = BlastParser()

        processed_gene_id = set()
        for hit in blast_parser.read_hit(table):
            genome_id, gene_id = hit.query_id.split('~')
            if gene_id in processed_gene_id:
                # Only consider the first hit as diamond/blast
                # tables are sorted by bitscore. In practice, few
                # genes will have multiple top hits.
                continue

            processed_gene_id.add(gene_id)
            scaffold_id = gene_id[0:gene_id.rfind('_')]

            subject_genome_id, subject_gene_id = hit.subject_id.split('~')
            self.profiles[genome_id].add_hit(gene_id,
                                             scaffold_id,
                                             subject_gene_id,
                                             subject_genome_id,
                                             taxonomy[subject_genome_id],
                                             hit.evalue,
                                             hit.perc_identity,
                                             hit.aln_length,
                                             hit.query_end - hit.query_start + 1)

    def write_genome_summary(self, output_file):
        """Summarize classification of each genome.

        Parameters
        ----------
        output_file : str
            Output file.
        """

        fout = open(output_file, 'w')
        fout.write('Genome id\t# scaffolds\t# genes\tCoding bases')
        for rank in Taxonomy.rank_labels:
            fout.write('\t' + rank + ': taxon')
            fout.write('\t' + rank + ': % of scaffolds')
            fout.write('\t' + rank + ': % of genes')
            fout.write('\t' + rank + ': % of coding bases')
            fout.write('\t' + rank + ': avg. e-value')
            fout.write('\t' + rank + ': avg. % identity')
            fout.write('\t' + rank + ': avg. align. length (aa)')
        fout.write('\n')

        sorted_genome_ids = alphanumeric_sort(self.profiles.keys())
        for genome_id in sorted_genome_ids:
            self.profiles[genome_id].write_genome_summary(fout)

        fout.close()
        
    def common_taxa(self, common_threshold, min_classified_per):
        """Get common taxa at each rank.
        
        Any rank in the dictionary that does not exist had 
        an insufficient number of classified genes to establish
        a set of common taxa.
        
        Parameters
        ----------
        common_threshold : float
            Percentage of classified genes for a taxon to be defined as common.
        min_classified_per : float
            Required percentage of classified genes to determine common taxa.

        Returns
        -------
        d[genome_id][rank] -> set of common taxa
            Support for taxon at each rank.
        """
        
        common_taxa = {}
        
        bin_report_dir = os.path.join(self.output_dir, 'bin_reports')
        for f in os.listdir(bin_report_dir):
            if not f.endswith('.gene.tsv'):
                continue
            
            genome_id = f[0:f.rfind('_genes.gene.tsv')]
            common_taxa[genome_id] = {}
        
            profile = self.read_scaffold_profile(genome_id, classified_genes=True)
            scaffold_stats = self.read_scaffold_stats(genome_id)
            
            # identify common taxa across scaffolds
            for rank in xrange(0, len(Taxonomy.rank_prefixes)):
                total_genes = 0
                genome_gene_count = defaultdict(int)
                for scaffold_id in profile:
                    _, _, _, num_genes, _ = scaffold_stats[scaffold_id]
                    total_genes += num_genes
                    
                    for taxon, stats in profile[scaffold_id][rank].iteritems():
                        _percent, gene_count, classified_genes = stats
                        genome_gene_count[taxon] += gene_count
                        
                total_classified_genes = sum(genome_gene_count.values())
                if total_genes == 0 or (total_classified_genes * 100.0 / total_genes < min_classified_per):
                    break
                    
                common_taxa[genome_id][rank] = set()
                for taxon, gene_count in genome_gene_count.iteritems():
                    if gene_count * 100.0 / total_classified_genes >= common_threshold:
                        common_taxa[genome_id][rank].add(taxon)
            
        return common_taxa
        
    def read_genome_profile(self):
        """Read taxonomic identification of each genome.

        Returns
        -------
        d[genome_id][rank] -> (taxon, percentage)
            Support for taxon at each rank.
        """
        
        profiles = {}
        
        genome_summary_file = os.path.join(self.output_dir, 'genome_summary.tsv')
        with open(genome_summary_file) as f:
            f.readline()
            
            for line in f:
                line_split = line.split('\t')
                genome_id = line_split[0]
                if genome_id.endswith('_genes'):
                    genome_id = genome_id[0:genome_id.rfind('_genes')]
                
                profiles[genome_id] = {}
                for r, i in enumerate(xrange(4, len(line_split), 7)):
                    taxon = line_split[i]
                    gene_support = float(line_split[i+2])
                    
                    if taxon == 'unclassified':
                        profiles[genome_id][r] = (Taxonomy.rank_prefixes[r], 0)
                    else:
                        profiles[genome_id][r] = (taxon, float(gene_support))

        return profiles
        
    def read_genome_taxonomy(self):
        """Read taxonomic identification of each genome.
        
        Taxonomy string has the form:
            d__<taxa> (<% genes>); p__<taxa> (<% genes>); ...; s__<taxa> (<% genes>)

        Returns
        -------
        d[genome_id] -> taxonomy info
            Taxonomic classification with percentage of supporting genes.
        """
        
        profiles = self.read_genome_profile()
        
        classification = {}
        for genome_id in profiles:
            taxa_str = []
            for r in xrange(0, len(Taxonomy.rank_prefixes)):
                taxa_str.append('%s (%.2f)' % profiles[genome_id].get(r, (Taxonomy.rank_prefixes[r], 0)))
                
            classification[genome_id] = ';'.join(taxa_str)

        return classification
        
    def read_scaffold_taxonomy(self):
        """Read taxonomic identification of each scaffold.
        
        Taxonomy string has the form:
            d__<taxa> (<% genes>); p__<taxa> (<% genes>); ...; s__<taxa> (<% genes>)

        Returns
        -------
        d[genome_id][scaffold_id] -> taxonomy info
            Taxonomic classification with percentage of supporting genes.
        """
        
        classification = defaultdict(dict)
        
        bin_report_dir = os.path.join(self.output_dir, 'bin_reports')
        for f in os.listdir(bin_report_dir):
            if not f.endswith('.scaffolds.tsv'):
                continue
                
            scaffold_summary_file = os.path.join(bin_report_dir, f)
            with open(scaffold_summary_file) as f:
                f.readline()
                
                for line in f:
                    line_split = line.split('\t')
                    scaffold_id = line_split[0]
                    genome_id = line_split[1]
                    if genome_id.endswith('_genes'):
                        genome_id = genome_id[0:genome_id.rfind('_genes')]
                    
                    taxa = []
                    for i in xrange(7, len(line_split), 5):
                        taxon = line_split[i]
                        gene_support = line_split[i+1]
                        
                        taxa.append('%s (%s)' % (taxon, gene_support))
                        
                    classification[genome_id][scaffold_id] = ';'.join(taxa)

        return classification

    def read_scaffold_profile(self, genome_id, classified_genes):
        """Read complete taxonomic profile.
        
        Parameters
        ----------
        genome_id : str
            Identifier of genome of interest.
        classified_genes : boolean
            Base profile percentages on just classified genes (True), or all genes (False).
            
        Returns
        -------
        d[scaffold_id][rank][taxon] -> (percentage, genes with classification, genes considered)
            Classification of scaffold across all ranks and taxa.
        """
        
        # read number of genes in each scaffold
        gene_count = {}
        with open(os.path.join(self.output_dir, 'bin_reports', genome_id + '_genes.scaffolds.tsv')) as fin:
            fin.readline()
            
            for line in fin:
                line_split = line.split('\t')
                
                scaffold_id = line_split[0]
                genes = int(line_split[5])
                
                gene_count[scaffold_id] = genes

        # read taxonomic assignment of each gene
        gene_taxonomy = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
        with open(os.path.join(self.output_dir, 'bin_reports', genome_id + '_genes.gene.tsv')) as fin:
            fin.readline()
            
            for line in fin:
                line_split = line.split('\t')
                
                gene_id = line_split[0]
                scaffold_id = gene_id[0:gene_id.rfind('_')]
                taxonomy = line_split[4].split(';')
                
                for r, t in enumerate(taxonomy):
                    if t == Taxonomy.rank_prefixes[r]:
                        continue
                        
                    gene_taxonomy[scaffold_id][r][t] += 1
                    
        # calculate percentages
        profile = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
        for scaffold_id in gene_taxonomy:
            for rank in xrange(0, len(Taxonomy.rank_prefixes)):
                total = gene_count[scaffold_id]
                if rank in gene_taxonomy[scaffold_id]:
                    if classified_genes:
                        total = sum(gene_taxonomy[scaffold_id][rank].values())
                    
                    for taxon, count in gene_taxonomy[scaffold_id][rank].iteritems():
                        profile[scaffold_id][rank][taxon] = (float(count) * 100.0 / total, count, total)
                else:
                    # no classification at this rank
                    if classified_genes:
                        total = 0
                    profile[scaffold_id][rank][Taxonomy.rank_prefixes[rank]] = (0.0, 0, total)
                    
        return profile
    
    def read_scaffold_stats(self, genome_id):
        """Read common statistics for scaffold.
        
        Parameters
        ----------
        genome_id : str
            Identifier of genome of interest.
            
        Returns
        -------
        d[scaffold_id] -> (length, GC, mean coverage, # genes, coding bases)
            Common statistics for each scaffold.
        """
        
        stats = {}
        with open(os.path.join(self.output_dir, 'bin_reports', genome_id + '_genes.scaffolds.tsv')) as fin:
            fin.readline()
            
            for line in fin:
                line_split = line.split('\t')
                
                scaffold_id = line_split[0]
                length = int(line_split[2])
                gc = float(line_split[3])
                mean_cov = float(line_split[4])
                genes = int(line_split[5])
                coding_bases = int(line_split[6])
                
                stats[scaffold_id] = (length, gc, mean_cov, genes, coding_bases)

        return stats
    
    def run(self, gene_files, 
                    stat_file, 
                    db_file, 
                    taxonomy_file, 
                    percent_to_classify, 
                    evalue, 
                    per_identity, 
                    per_aln_len,
                    tmpdir):
        """Create taxonomic profiles for a set of genomes.

        Parameters
        ----------
        gene_files : list of str
            Fasta files of called genes to process.
        stat_file : str
            File with statistics for individual scaffolds.
        db_file : str
            Database of reference genes.
        taxonomy_file : str
            File containing GreenGenes taxonomy strings for reference genomes.
        percent_to_classify : float
            Minimum percentage of genes to assign scaffold to a taxon [0, 100].
        evalue : float
            E-value threshold used to identify homologs.
        per_identity: float
            Percent identity threshold used to identify homologs [0, 100].
        per_aln_len : float
            Percent coverage of query sequence used to identify homologs [0, 100].
        tmpdir : str
            Directory to use for temporary files.
        """

        # read statistics file
        self.logger.info('Reading scaffold statistics.')
        scaffold_stats = ScaffoldStats()
        scaffold_stats.read(stat_file)

        # concatenate gene files
        self.logger.info('Appending genome identifiers to all gene identifiers.')
        diamond_output_dir = os.path.join(self.output_dir, 'diamond')
        make_sure_path_exists(diamond_output_dir)

        gene_file = os.path.join(diamond_output_dir, 'genes.faa')
        concatenate_gene_files(gene_files, gene_file)

        # read taxonomy file
        self.logger.info('Reading taxonomic assignment of reference genomes.')

        t = Taxonomy()
        taxonomy = t.read(taxonomy_file)
        if not t.validate(taxonomy, 
                            check_prefixes=True, 
                            check_ranks=True, 
                            check_hierarchy=False, 
                            check_species=False,
                            check_group_names=False,
                            check_duplicate_names=False, 
                            report_errors=True):
            self.logger.error('Invalid taxonomy file.')
            sys.exit()
            
        # record length and number of genes in each scaffold
        for aa_file in gene_files:
            genome_id = remove_extension(aa_file)
            self.profiles[genome_id] = Profile(genome_id, percent_to_classify, taxonomy)

            for seq_id, seq in seq_io.read_seq(aa_file):
                scaffold_id = seq_id[0:seq_id.rfind('_')]
                self.profiles[genome_id].genes_in_scaffold[scaffold_id] += 1
                self.profiles[genome_id].coding_bases[scaffold_id] += len(seq) * 3  # length in nucleotide space

        # run diamond and create taxonomic profile for each genome
        self.logger.info('Running diamond blastp with %d processes (be patient!)' % self.cpus)

        diamond = Diamond(self.cpus)
        diamond_table_out = os.path.join(diamond_output_dir, 'diamond_hits.tsv')
        diamond.blastp(gene_file, 
                        db_file, 
                        evalue, 
                        per_identity, 
                        per_aln_len, 
                        1, 
                        False,
                        diamond_table_out, 
                        output_fmt='standard', 
                        tmp_dir=tmpdir)
               
        # create taxonomic profile for each genome
        self.logger.info('Creating taxonomic profile for each genome.')
        self.taxonomic_profiles(diamond_table_out, taxonomy)

        # write out taxonomic profile
        self.logger.info('Writing taxonomic profile for each genome.')
        report_dir = os.path.join(self.output_dir, 'bin_reports')
        make_sure_path_exists(report_dir)
        for aa_file in gene_files:
            genome_id = remove_extension(aa_file)
            profile = self.profiles[genome_id]

            scaffold_summary_out = os.path.join(report_dir, genome_id + '.scaffolds.tsv')
            profile.write_scaffold_summary(scaffold_stats, scaffold_summary_out)

            gene_summary_out = os.path.join(report_dir, genome_id + '.gene.tsv')
            profile.write_gene_summary(gene_summary_out, seq_io.read(aa_file))

            genome_profile_out = os.path.join(report_dir, genome_id + '.profile.tsv')
            profile.write_genome_profile(genome_profile_out)

        # create summary report for all genomes
        genome_summary_out = os.path.join(self.output_dir, 'genome_summary.tsv')
        self.write_genome_summary(genome_summary_out)

        # create Krona plot based on classification of scaffolds
        self.logger.info('Creating Krona plot for each genome.')
        krona_profiles = defaultdict(lambda: defaultdict(int))
        for genome_id, profile in self.profiles.iteritems():
            seq_assignments = profile.classify_seqs()

            for seq_id, classification in seq_assignments.iteritems():
                taxa = []
                for r in xrange(0, len(Taxonomy.rank_labels)):
                    taxa.append(classification[r][0])

                krona_profiles[genome_id][';'.join(taxa)] += profile.genes_in_scaffold[seq_id]

        krona = Krona()
        krona_output_file = os.path.join(self.output_dir, 'gene_profiles.scaffolds.html')
        krona.create(krona_profiles, krona_output_file)

        # create Krona plot based on best hit of each gene
        krona_profiles = defaultdict(lambda: defaultdict(int))

        for aa_file in gene_files:
            genome_id = remove_extension(aa_file)

            profile = self.profiles[genome_id]
            for gene_id, _seq in seq_io.read_seq(aa_file):

                taxa_str = Taxonomy.unclassified_taxon
                if gene_id in profile.gene_hits:
                    taxa_str, _hit_info = profile.gene_hits[gene_id]

                krona_profiles[genome_id][taxa_str] += 1

        krona_output_file = os.path.join(self.output_dir, 'gene_profiles.genes.html')
        krona.create(krona_profiles, krona_output_file)
        
    def filter(self, 
                consensus_taxon_threshold, 
                trusted_scaffold_threshold, 
                common_taxa_threshold, 
                congruent_scaffold_threshold,
                min_classified_per_threshold,
                min_classified_threshold,
                consensus_scaffold_threshold,
                output_file):
        """Filter scaffolds with divergent taxonomic classification.
        
        Parameters
        ----------
        consensus_taxon_threshold : float
            Threshold for accepting a consensus taxon.
        trusted_scaffold_threshold : float
            Threshold for treating a scaffold as trusted.
        common_taxa_threshold : float
            Threshold for treating a taxon as common.
        congruent_scaffold_threshold : float
            Threshold for treating a scaffold as congruent.
        min_classified_per_threshold : float
            Minimum percentage of genes with a classification to filter a scaffold.
        min_classified_threshold : int
            Minimum number of classified genes required to filter a scaffold.
        consensus_scaffold_threshold : float
            Threshold of consensus taxon for filtering a scaffold.
        output_file : str
            File to write filtered scaffolds.
        """
        
        # filter scaffolds with divergent taxonomic classifications
        self.logger.info('Identifying scaffolds with divergent taxonomic classifications.')
        fout = open(output_file, 'w')
        fout.write('# Taxon filtering with RefineM v%s\n' % version())
        fout.write('# consensus_taxon_threshold: %.2f\n' % consensus_taxon_threshold)
        fout.write('# trusted_scaffold_threshold: %.2f\n' % trusted_scaffold_threshold)
        fout.write('# common_taxa_threshold: %.2f\n' % common_taxa_threshold)
        fout.write('# congruent_scaffold_threshold: %.2f\n' % congruent_scaffold_threshold)
        fout.write('# min_classified_per_threshold: %.2f\n' % min_classified_per_threshold)
        fout.write('# min_classified_threshold: %.2f\n' % min_classified_threshold)
        fout.write('Scaffold id\tGenome id\t# classified scaffolds')
        fout.write('\tConsensus taxon\tGenome support\tScaffold support')
        fout.write('\t# trusted scaffolds\tCommon taxa\tCommon taxa support')
        fout.write('\tScaffold taxon\tScaffold support')
        fout.write('\tLength (bp)\t# genes\t# classified genes\tGC\tMean coverage\n')
        
        bin_report_dir = os.path.join(self.output_dir, 'bin_reports')
        for f in os.listdir(bin_report_dir):
            if not f.endswith('.gene.tsv'):
                continue
            
            genome_id = f[0:f.rfind('_genes.gene.tsv')]
            profile = self.read_scaffold_profile(genome_id, classified_genes=True)
            scaffold_stats = self.read_scaffold_stats(genome_id)
            
            for rank in xrange(0, len(Taxonomy.rank_prefixes)):
                # determine consensus taxon for genome
                genome_consensus = defaultdict(int)
                genome_classified_genes = 0
                for scaffold_id in profile:
                    for taxon, stats in profile[scaffold_id][rank].iteritems():
                        _, gene_count, classified_genes = stats
                        
                        genome_consensus[taxon] += gene_count
                        
                    genome_classified_genes += classified_genes
                    
                if genome_classified_genes == 0:
                    # no classifed genes at this rank
                    break
                    
                consensus_taxon, consensus_genes = sorted(genome_consensus.items(), key=operator.itemgetter(1), reverse=True)[0]
                consensus_support = float(consensus_genes) * 100.0 / genome_classified_genes
                if consensus_support < consensus_taxon_threshold:
                    # stop filtering
                    break
                    
                # identify trusted scaffolds based on consensus taxon
                trusted_scaffolds = set()
                trusted_gene_profile = defaultdict(float)
                trusted_classified_genes = 0
                for scaffold_id in profile:
                    support, _, _ = profile[scaffold_id][rank].get(consensus_taxon, [0, 0, 0])
                    if support > trusted_scaffold_threshold:
                        trusted_scaffolds.add(scaffold_id)
                        
                        for taxon, stats in profile[scaffold_id][rank].iteritems():
                            _, gene_count, classified_genes = stats
                            trusted_gene_profile[taxon] += gene_count

                        trusted_classified_genes += classified_genes

                # identify common taxa across trusted scaffolds
                common_taxa = set()
                for taxon, gene_count in trusted_gene_profile.iteritems():
                    if gene_count * 100.0 / trusted_classified_genes > common_taxa_threshold:
                        common_taxa.add(taxon)
     
                # filter scaffolds that are incongruent with list of common taxa
                filtered_scaffolds = []
                for scaffold_id in profile:
                    length, gc, mean_cov, num_genes, _coding_bases = scaffold_stats[scaffold_id]
                    scaffold_consensus_support, _, _ = profile[scaffold_id][rank].get(consensus_taxon, [0, 0, 0])
                    classified_genes = profile[scaffold_id][rank].values()[0][2]

                    if classified_genes < max(min_classified_threshold, num_genes * min_classified_per_threshold/100):
                        continue
                
                    congruent_gene_count = 0
                    scaffold_taxon = consensus_taxon
                    scaffold_taxon_support = scaffold_consensus_support
                    for taxon, stats in profile[scaffold_id][rank].iteritems():
                        support, gene_count, classified_genes = stats
                        
                        if taxon in common_taxa:                            
                            congruent_gene_count += gene_count
                            
                        if support > scaffold_taxon_support:
                            scaffold_taxon_support = support
                            scaffold_taxon = taxon
                            
                    if classified_genes == 0:
                        # can not filter a scaffold with no classified genes
                        # at the current rank
                        continue
                            
                    congruent_per = congruent_gene_count * 100.0 / classified_genes
                    if congruent_per <= congruent_scaffold_threshold and scaffold_taxon_support > consensus_scaffold_threshold:
                        common_taxa_str = ','.join(sorted(list(common_taxa)))
                        
                        fout.write('%s\t%s\t%d' % (scaffold_id, genome_id, len(profile)))
                        fout.write('\t%s\t%.2f\t%.2f' % (consensus_taxon, consensus_support, scaffold_consensus_support))
                        fout.write('\t%d\t%s\t%.2f' % (len(trusted_scaffolds), common_taxa_str, congruent_per))
                        fout.write('\t%s\t%.2f' % (scaffold_taxon, scaffold_taxon_support))
                        fout.write('\t%d\t%d\t%d\t%.2f\t%.2f\n' % (length, num_genes, classified_genes, gc, mean_cov))
                        
                        filtered_scaffolds.append(scaffold_id)
                        
                # remove filtered scaffolds before considering next rank
                for scaffold_id in filtered_scaffolds:
                    profile.pop(scaffold_id)

        fout.close()
   
    def filter_deprecated(self, genome_threshold, min_scaffold_agreement, max_scaffold_disagreement, min_classified_per, output_file):
        """Filter scaffolds with divergent taxonomic classification.
        
        Parameters
        ----------
        genome_threshold : float
            Threshold for accepting taxonomic classification of genome.
        min_scaffold_agreement : float
            Minimum percentage of genes congruent with genome classification to retain scaffold.
        max_scaffold_disagreement : float
            Maximum percentage of genes supporting an alternative taxon to retain scaffold.
        min_classified_per : float
            Minimum percentage of genes with a classification to filter a scaffold.
        output_file : str
            File to write filtered scaffolds.
        """
        
        # read taxonomic profiles for genomes: d[genome_id][rank] -> (taxon, support)
        self.logger.info('Reading genome profiles.')
        genome_profiles = self.read_genome_profile()
        
        # filter scaffolds with divergent taxonomic classifications
        self.logger.info('Identifying scaffolds with divergent taxonomic classifications.')
        fout = open(output_file, 'w')
        fout.write('# Taxon filtering with RefineM v%s\n' % version())
        fout.write('# min_scaffold_agreement: %.2f\n' % min_scaffold_agreement)
        fout.write('# max_scaffold_disagreement: %.2f\n' % max_scaffold_disagreement)
        fout.write('# min_classified_per: %.2f\n' % min_classified_per)
        fout.write('Scaffold id\tGenome id\tGenome taxon\tGenome support\tScaffold support\tScaffold taxon\tScaffold support\tLength (bp)\t# genes\t# classified genes\tGC\tMean coverage\n')
        
        bin_report_dir = os.path.join(self.output_dir, 'bin_reports')
        for f in os.listdir(bin_report_dir):
            if not f.endswith('.gene.tsv'):
                continue
            
            genome_id = f[0:f.rfind('_genes.gene.tsv')]
            profile = self.read_scaffold_profile(genome_id, classified_genes=True)
            
            scaffold_stats = self.read_scaffold_stats(genome_id)

            for scaffold_id in profile:
                length, gc, mean_cov, num_genes, _coding_bases = scaffold_stats[scaffold_id]
                
                for rank in profile[scaffold_id]:
                    genome_taxon, genome_support = genome_profiles[genome_id][rank]
                    
                    if genome_support >= genome_threshold:
                        scaffold_support, _, _ = profile[scaffold_id][rank].get(genome_taxon, [0, 0, 0])
                        
                        # determine classifcation of scaffold with the most support
                        scaffold_taxon_support = scaffold_support
                        scaffold_taxon = genome_taxon
                        for taxon, stats in profile[scaffold_id][rank].iteritems():
                            support, _, classified_genes = stats
                            if support > scaffold_taxon_support:
                                scaffold_taxon = taxon
                                scaffold_taxon_support = support

                        if (float(classified_genes) * 100 / num_genes) < min_classified_per:
                            # insufficient number of classified genes to filter 
                            # scaffold based on taxonomic classifcation
                            break
   
                        if (scaffold_support < min_scaffold_agreement 
                            or (scaffold_taxon_support > max_scaffold_disagreement and scaffold_taxon != genome_taxon)):
                            fout.write('%s\t%s' % (scaffold_id, genome_id))
                            fout.write('\t%s\t%.2f\t%.2f' % (genome_taxon, genome_support, scaffold_support))
                            fout.write('\t%s\t%.2f' % (scaffold_taxon, scaffold_taxon_support))
                            fout.write('\t%d\t%d\t%d\t%.2f\t%.2f\n' % (length, num_genes, classified_genes, gc, mean_cov))
                            
                            break
                    else:
                        # stop considering taxa once we reach a rank without sufficient support
                        break
        fout.close()

class Profile(object):
    """Profile of hits to reference genomes."""

    def __init__(self, genome_id, percent_to_classify, taxonomy):
        """Initialization.

        Parameters
        ----------
        genome_id : str
            Unique identify of genome.
        percent_to_classify : float
            Minimum percentage of genes to assign scaffold to a taxon [0, 100].
        taxonomy : d[ref_genome_id] -> [domain, phylum, ..., species]
            Taxonomic assignment of each reference genome.
        """

        self.percent_to_classify = percent_to_classify / 100.0

        self.unclassified = Taxonomy.unclassified_rank

        self.genome_id = genome_id
        self.taxonomy = taxonomy

        self.TaxaInfo = namedtuple('TaxaInfo', """evalue
                                                perc_identity
                                                aln_length
                                                num_seqs
                                                num_genes
                                                num_basepairs""")

        # track hits at each rank: dict[scaffold_id][rank][taxa] -> [HitInfo, ...]
        self.HitInfo = namedtuple('HitInfo', """subject_genome_id
                                                subject_gene_id
                                                evalue
                                                perc_identity
                                                aln_length
                                                query_aln_length""")
        self.hits = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        # hit information for individual genes: d[gene_id] -> [taxonomy_str, HitInfo]
        self.gene_hits = {}

        # number of coding bases in scaffold in nucleotide space
        self.coding_bases = defaultdict(int)

        # number of genes in each scaffold
        self.genes_in_scaffold = defaultdict(int)

    def add_hit(self,
                query_gene_id, query_scaffold_id,
                subject_gene_id, subject_genome_id,
                tax_list, evalue, perc_identity,
                aln_length, query_aln_length):
        """Add hit to profile.

        Parameters
        ----------
        query_gene_id : str
            Unique identifier of query gene.
        query_scaffold_id : str
            Unique identifier of query scaffold.
        subject_gene_id : str
            Unique identifier of subject gene.
        subject_genome_id : str
            Unique identifier of subject gene.
        tax_list : list
            List indication taxa at each taxonomic rank.
        evalue : float
            E-value of hit.
        per_identity: float
            Percent identity of hit.
        aln_length: int
            Alignment length of hit.
        query_aln_length : int
            Length of query sequence in alignment.
        """

        hit_info = self.HitInfo(subject_genome_id, subject_gene_id,
                                evalue, perc_identity,
                                aln_length, query_aln_length)

        self.gene_hits[query_gene_id] = [';'.join(tax_list), hit_info]

        d = self.hits[query_scaffold_id]
        for i, taxa in enumerate(tax_list):
            d[i][taxa].append(hit_info)

    def classify_seqs(self):
        """Classify scaffold.

        Scaffold are classified using a majority vote
        over all genes with a valid hit. If less than 20% of
        genes have a valid hit, the scaffold is considered unclassified.
        Classification is performed from the highest (domain)
        to lowest (species) rank. If a rank is taxonomically
        inconsistent with a higher ranks classification, this
        rank and all lower ranks are set to unclassified.

        Returns
        -------
        dict : d[scaffold_id][rank] -> [taxon, HitInfo]
            Classification of each scaffold along with summary statistics
            of hits to the specified taxon.
        """

        expected_parent = Taxonomy().taxonomic_consistency(self.taxonomy)

        # classify each scaffold using a majority vote
        seq_assignments = defaultdict(lambda: defaultdict(list))
        for seq_id, rank_hits in self.hits.iteritems():
            parent_taxa = None
            for rank in xrange(0, len(Taxonomy.rank_prefixes)):
                taxa = max(rank_hits[rank], key=lambda x: len(rank_hits[rank][x]))
                count = len(rank_hits[rank][taxa])
                
                if (taxa != Taxonomy.rank_prefixes[rank]
                    and (count >= self.percent_to_classify * self.genes_in_scaffold[seq_id])
                    and (rank == 0 or expected_parent[taxa] == parent_taxa)):
                        seq_assignments[seq_id][rank] = [taxa, rank_hits[rank][taxa]]
                        parent_taxa = taxa
                else:
                    # set to unclassified at all lower ranks
                    for r in xrange(rank, len(Taxonomy.rank_prefixes)):
                        seq_assignments[seq_id][r] = [self.unclassified, None]
                    break

        # identify scaffold with no hits
        for seq_id in self.genes_in_scaffold:
            if seq_id not in seq_assignments:
                for rank in xrange(0, len(Taxonomy.rank_prefixes)):
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
            Relative abundance of taxa at a given rank determined
            from the classification of each scaffolds and weighted by
            the number of genes in each scaffold.
        dict : d[rank][taxa] -> TaxaInfo
           Statistics for each taxa.
        """

        seq_assignments = self.classify_seqs()

        total_genes = sum(self.genes_in_scaffold.values())

        profile = defaultdict(lambda: defaultdict(float))
        stats = defaultdict(dict)

        for r in xrange(0, len(Taxonomy.rank_labels)):
            num_seqs = defaultdict(int)
            num_genes = defaultdict(int)
            num_basepairs = defaultdict(int)
            hit_stats = defaultdict(list)
            for seq_id, data in seq_assignments.iteritems():
                taxa, hit_info = data[r]
                profile[r][taxa] += float(self.genes_in_scaffold[seq_id]) / total_genes
                num_seqs[taxa] += 1
                num_genes[taxa] += self.genes_in_scaffold[seq_id]
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
                                            num_genes[taxa],
                                            num_basepairs[taxa])

            stats[r][self.unclassified] = self.TaxaInfo(None,
                                                     None,
                                                     None,
                                                     num_seqs[self.unclassified],
                                                     num_genes[self.unclassified],
                                                     num_basepairs[self.unclassified])

        return profile, stats

    def write_genome_summary(self, fout):
        """Write profile of most abundant taxon at each rank.

        Parameters
        ----------
        fout : output stream
            Output stream.
        """

        profile, stats = self.profile()

        total_seqs = len(self.genes_in_scaffold)
        total_genes = sum(self.genes_in_scaffold.values())
        total_coding_bases = sum(self.coding_bases.values())

        fout.write('%s\t%d\t%d\t%d' % (self.genome_id, total_seqs, total_genes, total_coding_bases))
        for r in xrange(0, len(Taxonomy.rank_labels)):
            if len(profile[r]) == 0:
                continue
                
            taxa, _percent = max(profile[r].iteritems(), key=operator.itemgetter(1))

            if taxa != self.unclassified:
                fout.write('\t%s\t%.2f\t%.2f\t%.2f\t%.2g\t%.2f\t%.2f' % (taxa,
                                                              stats[r][taxa].num_seqs * 100.0 / total_seqs,
                                                              stats[r][taxa].num_genes * 100.0 / total_genes,
                                                              stats[r][taxa].num_basepairs * 100.0 / total_coding_bases,
                                                              stats[r][taxa].evalue,
                                                              stats[r][taxa].perc_identity,
                                                              stats[r][taxa].aln_length))
            else:
                fout.write('\t%s\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s' % (taxa,
                                                          stats[r][taxa].num_seqs * 100.0 / total_seqs,
                                                          stats[r][taxa].num_genes * 100.0 / total_genes,
                                                          stats[r][taxa].num_basepairs * 100.0 / total_coding_bases,
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
        for rank in Taxonomy.rank_labels:
            if rank != Taxonomy.rank_labels[0]:
                fout.write('\t')

            fout.write(rank + ': taxon')
            fout.write('\t' + rank + ': % scaffolds')
            fout.write('\t' + rank + ': % genes')
            fout.write('\t' + rank + ': % coding bps')
            fout.write('\t' + rank + ': avg. e-value')
            fout.write('\t' + rank + ': avg. % identity')
            fout.write('\t' + rank + ': avg. align. length (aa)')
        fout.write('\n')

        # sort profiles in descending order of abundance
        profile, stats = self.profile()

        sorted_profiles = {}
        max_taxa = 0
        for r in xrange(0, len(Taxonomy.rank_labels)):
            sorted_profile = sorted(profile[r].items(), key=operator.itemgetter(1))
            sorted_profile.reverse()

            sorted_profiles[r] = sorted_profile

            if len(sorted_profiles) > max_taxa:
                max_taxa = len(sorted_profiles)

        # write out table
        total_seqs = len(self.genes_in_scaffold)
        total_genes = sum(self.genes_in_scaffold.values())
        total_coding_bases = sum(self.coding_bases.values())

        for i in xrange(0, max_taxa):
            for r in xrange(0, len(Taxonomy.rank_labels)):
                if r != 0:
                    fout.write('\t')

                if len(sorted_profiles[r]) > i:
                    taxa, _percent = sorted_profiles[r][i]

                    if taxa != self.unclassified:
                        fout.write('%s\t%.2f\t%.2f\t%.2f\t%.2g\t%.2f\t%.2f' % (taxa,
                                                                      stats[r][taxa].num_seqs * 100.0 / total_seqs,
                                                                      stats[r][taxa].num_genes * 100.0 / total_genes,
                                                                      stats[r][taxa].num_basepairs * 100.0 / total_coding_bases,
                                                                      stats[r][taxa].evalue,
                                                                      stats[r][taxa].perc_identity,
                                                                      stats[r][taxa].aln_length))
                    else:
                        fout.write('%s\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s' % (taxa,
                                                                  stats[r][taxa].num_seqs * 100.0 / total_seqs,
                                                                  stats[r][taxa].num_genes * 100.0 / total_genes,
                                                                  stats[r][taxa].num_basepairs * 100.0 / total_coding_bases,
                                                                  'na',
                                                                  'na',
                                                                  'na'))
                else:
                    fout.write('\t\t\t\t\t')

            fout.write('\n')

    def write_scaffold_summary(self, scaffold_stats, output_file):
        """Summarize classification of each scaffold.

        Parameters
        ----------
        scaffold_stats : ScaffoldStats
            Class containing statistics for scaffolds.
        output_file : str
            Output file.
        """

        seq_assignments = self.classify_seqs()

        fout = open(output_file, 'w')
        fout.write('Scaffold id')
        fout.write('\tGenome id\tLength (bp)\tGC\tMean coverage')
        fout.write('\t# genes\tCoding bases (nt)')
        for rank in Taxonomy.rank_labels:
            fout.write('\t' + rank + ': taxa')
            fout.write('\t' + rank + ': % genes')
            fout.write('\t' + rank + ': avg. e-value')
            fout.write('\t' + rank + ': avg. % identity')
            fout.write('\t' + rank + ': avg. align. length (aa)')
        fout.write('\n')

        for seq_id in seq_assignments:
            fout.write('%s\t%s\t%.2f\t%d\t%d' % (seq_id,
                                       scaffold_stats.print_stats(seq_id),
                                       mean(scaffold_stats.coverage(seq_id)),
                                       self.genes_in_scaffold[seq_id],
                                       self.coding_bases[seq_id]))

            for r in xrange(0, len(Taxonomy.rank_labels)):
                taxa, hit_info = seq_assignments[seq_id][r]

                if taxa != self.unclassified:
                    avg_evalue = mean([x.evalue for x in hit_info])
                    avg_perc_identity = mean([x.perc_identity for x in hit_info])
                    avg_aln_length = mean([x.aln_length for x in hit_info])

                    hit_str = '%.2f' % (len(hit_info) * 100.0 / self.genes_in_scaffold[seq_id])
                    fout.write('\t%s\t%s\t%.2g\t%.2f\t%.2f' % (taxa,
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

    def write_gene_summary(self, output_file, gene_seqs):
        """Summarize classification of each gene.

        Parameters
        ----------
        output_file : str
            Output file.
        gene_seqs : d[gene_id] -> amino acid sequence
            Amino acid sequence of each gene.
        """

        fout = open(output_file, 'w')
        fout.write('Gene id\tCoding bases (nt)\tSubject genome id\tSubject gene id\tTaxonomy\te-value\t% identity\talign. length (aa)\t% query aligned\tQuery sequence\n')

        for gene_id, data in self.gene_hits.iteritems():
            taxonomy, hit_info = data

            seq = gene_seqs[gene_id]
            if seq[-1] == '*':
                seq = seq[0:-1]

            fout.write('%s\t%d\t%s\t%s\t%s\t%.2g\t%.2f\t%d\t%.2f\t%s\n' % (gene_id,
                                                                     len(seq) * 3,
                                                                     hit_info.subject_genome_id,
                                                                     hit_info.subject_gene_id,
                                                                     taxonomy,
                                                                     hit_info.evalue,
                                                                     hit_info.perc_identity,
                                                                     hit_info.aln_length,
                                                                     hit_info.query_aln_length * 100.0 / len(seq),
                                                                     seq))
