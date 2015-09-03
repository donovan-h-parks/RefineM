#!/usr/bin/env python

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

__prog_name__ = 'make_database'
__prog_desc__ = 'create dereplicated database of genes from reference genomes'

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os
import sys
import argparse
import datetime
import string
import random
import tempfile
import shutil
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import make_sure_path_exists, remove_extension
from biolib.taxonomy import Taxonomy
from biolib.external.execute import check_dependencies
from biolib.misc.time_keeper import TimeKeeper


class MakeDatabase(object):
    """Make a dereplicated database of genes.

    Dereplication is done between genes within a named taxonomic
    group (e.g., genomes in the same genus) and is based on the
    average amino acid identity (AAI) between genes. Groups with large
    numbers of taxa can take an excessive amount of time to
    process so are subsampled to a specific number of taxa.
    Subsampling is done in a manor which aims to retain
    phylogenetic diversity and thus helps ensures a good
    distribution of genes within the group. Care is taken
    to ensure type strains are retained during dereplication.

    Note: this script is tailored to IMG in that it assumes
    a certain directory structure and file extensions. It also
    corrects a number of common issues with IMG genomes:
      - non-ascii characters in fasta header lines
      - hyphens at the start of some protein sequences
    """

    def __init__(self):
        """Initialize."""

        check_dependencies(['comparem', 'diamond', 'makeblastdb'])

        self.underclassified = 'underclassified'

        self.rank_prefixes = Taxonomy.rank_prefixes
        self.rank_index = Taxonomy.rank_index
        self.rank_labels = Taxonomy.rank_labels

        self.time_keeper = TimeKeeper()

    def read_taxonomy(self, input_taxonomy):
        """Read taxonomy file.

        Taxonomy file should have the following format:
            <genome_id>\t<taxonomy_str>

            where taxonomy_str is in GreenGenes format:
                d__Bacteria;p__Firmicutes;...

        Parameters
        ----------
        input_taxonomy : str
            Taxonomy file.

        Returns
        -------
        dict
            Taxonomy for each genome id.
        """

        taxonomy = {}
        for line in open(input_taxonomy):
            line_split = line.split('\t')

            taxonomy[line_split[0]] = [x.strip() for x in line_split[1].split(';')]

        return taxonomy

    def read_type_strain(self, type_strain_file):
        """Read type strain file.

        The type strain file should have the following format:
            <genome_id>\t<genome_name>

        Parameters
        ----------
        type_strain_file : str
            File specifying type strains.

        Returns
        -------
        set
            Set of all genome ids specified as type strains.
        """

        type_strains = set()
        for line in open(type_strain_file):
            line_split = line.split('\t')
            type_strains.add(line_split[0])

        return type_strains

    def select_taxa(self, genome_list, taxonomy, type_strains, max_taxa):
        """Select subset of genomes with a good distribution across named groups.

        Groups genomes into named groups and subsamples evenly across
        these groups. Ideally, genomes would be grouped into species, but
        some genomes may not have a species identifier. Such genomes are
        assigned to the most specific named group possible. Any genome
        marked as a type strain will be retained.

        Parameters
        ----------
        genome_list : iterable of genome ids
            Genomes to subsample.
        taxonomy : d[genome_id] -> [domain, ..., species]
            Taxonomy of each genome.
        type_strains : iterable
            Genome identifiers of type strains.
        max_taxa : int
            Number of genomes to retain.

        Returns
        -------
        iterable
            Subsampled list of genomes.
        """

        if len(genome_list) <= max_taxa:
            return genome_list

        reduced_genome_list = []

        # group genomes into the most specific named groups possible
        groups = defaultdict(set)
        for genome_id in genome_list:
            # add in type strains regardless of taxonomy
            if genome_id in type_strains:
                reduced_genome_list.append(genome_id)
                continue

            # get first classified rank
            for rank_index in xrange(self.rank_index['s__'], -1, -1):
                taxa = taxonomy[genome_id][rank_index]
                if taxa != self.rank_prefixes[rank_index]:
                    break

            groups[taxa].add(genome_id)

        # sample genomes from each named group
        while len(reduced_genome_list) < max_taxa:
            genomes_to_select = max_taxa - len(reduced_genome_list)
            genomes_per_group = max(genomes_to_select / len(groups), 1)
            for taxa, genome_ids in groups.iteritems():
                selected_genomes = random.sample(genome_ids, min(len(genome_ids), genomes_per_group))
                groups[taxa] = genome_ids.difference(selected_genomes)

                reduced_genome_list.extend(selected_genomes)

                if len(reduced_genome_list) == max_taxa:
                    break  # special case where we are adding single genomes from each group

        return reduced_genome_list

    def write_gene_file(self, gene_out, gene_dir, genome_list, taxonomy, genes_to_ignore):
        """Write genes to output stream.

        Parameters
        ----------
        gene_out : stream
            Output stream.
        gene_dir : str
            Directory containing called genes in amino acid space.
        genome_list : iterable
            Genomes to process.
        genes_to_ignore : set
            Genes which should not be written to file.
        """

        genes_kept = 0
        for genome_id in genome_list:
            genome_gene_file = os.path.join(gene_dir, genome_id + '.faa')
            if not os.path.exists(genome_gene_file):
                print '[WARNING] Missing gene file for genome %s.' % genome_gene_file
                continue

            if os.stat(genome_gene_file).st_size == 0:
                print '[WARNING] Gene file is empty for genome %s.' % genome_gene_file
                continue

            for gene_id, seq, annotation in seq_io.read_fasta_seq(genome_gene_file, keep_annotation=True):
                if gene_id in genes_to_ignore:
                    continue

                # IMG headers sometimes contain non-ascii characters which cause
                # problems with BLAST and DIAMOND so there are explicitly filtered out
                annotation = filter(lambda x: x in string.printable, annotation)

                # a few IMG genomes contain protein sequences which start with a hyphen
                if seq[0] == '-':
                    seq = seq[1:]

                gene_out.write('>' + gene_id + ' ' + annotation + '\n')
                gene_out.write(seq + '\n')
                genes_kept += 1

        return genes_kept

    def img_gene_id_to_scaffold_id(self, genome_dir, genome_id, output_dir):
        """Modify IMG gene ids to format which explicitly gives scaffold names.

        For downstream processing it is often necessary to know which scaffold
        a gene is contained on. IMG uses unique identifiers for genes. As such,
        these are changed to the following format:

        <scaffold_id>_<gene #> <annotation> [IMG gene id]

        Parameters
        ----------
        genome_dir : str
            Directory with files for genome.
        genome_id : str
            Unique identifier of genome.
        output_dir : float
            Directory to contain modified fasta files.
        """

        # determine source scaffold for each gene
        gene_id_to_scaffold_id = {}
        gene_number = defaultdict(int)
        for line in open(os.path.join(genome_dir, genome_id + '.gff')):
            if line[0] == '#':
                continue

            line_split = line.split('\t')
            scaffold_id = line_split[0]
            info = line_split[8]
            if info != '':  # this will be empty for non-protein coding genes
                gene_id = info.split(';')[0].replace('ID=', '')

                gene_number[scaffold_id] += 1
                gene_id_to_scaffold_id[gene_id] = scaffold_id + '_' + str(gene_number[scaffold_id])

        # write out gene file with modified identifiers
        genome_gene_file = os.path.abspath(os.path.join(genome_dir, genome_id + '.genes.faa'))

        fout = open(os.path.join(output_dir, genome_id + '.faa'), 'w')
        for gene_id, seq, annotation in seq_io.read_fasta_seq(genome_gene_file, keep_annotation=True):

            annotation = annotation[annotation.find(' ') + 1:]  # remove additional gene id from annotation
            annotation += ' [IMG Gene ID: ' + gene_id + ']'  # append IMG gene id for future reference

            fout.write('>' + gene_id_to_scaffold_id[gene_id] + ' ' + annotation + '\n')
            fout.write(seq + '\n')
        fout.close()

    def amend_gene_identifies(self, gene_dir, output_dir):
        """Modify gene ids to include source genome id.

        The following format is used:
          <gene_id>~<genome_id>

        Parameters
        ----------
        gene_dir : str
            Directory with fasta files containing protein sequences.
        output_dir : float
            Directory to contain modified fasta files.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for f in os.listdir(gene_dir):
            gf = os.path.join(gene_dir, f)
            genome_id = remove_extension(gf)

            aa_file = os.path.join(output_dir, genome_id + '.faa')
            fout = open(aa_file, 'w')
            for seq_id, seq, annotation in seq_io.read_fasta_seq(gf, keep_annotation=True):
                fout.write('>' + seq_id + '~' + genome_id + ' ' + annotation + '\n')
                if seq[-1] == '*':
                    seq = seq[0:-1]
                fout.write(seq + '\n')
            fout.close()

    def filter_aai(self, tmp_dir, gene_dir, ammended_gene_dir, per_identity, per_aln_len, cpus):
        """Filter genes with similar amino acid identity.

        Parameters
        ----------
        tmp_dir : str
            Temporary directory for storing results.
        gene_dir : str
            Directory with fasta files containing protein sequences.
        ammended_gene_dir : str
            Directory to store protein sequences with ammended gene ids.
        per_identity : float
            Percent identity for subsampling similar genes.
        per_aln_len : float
            Percent alignment length for subsampling similar genes.
        cpus : int
            Number of cpus to use.

        Returns
        -------
        genes_to_remove : set
            Unique identifiers of genes to filter.
        """

        rblast_dir = os.path.join(tmp_dir, 'rblast')
        os.system('comparem rblast -e 1e-10 -p %d -c %d %s %s' % (per_identity, cpus, gene_dir, rblast_dir))
        aai_dir = os.path.join(tmp_dir, 'aai')
        os.system('comparem aai -p %d -a %d -c %d %s %s' % (per_identity, per_aln_len, cpus, rblast_dir, aai_dir))

        # identify homologs to be filtered
        print ''
        print '  Identifying homologs to be filtered.'
        shared_genes_dir = os.path.join(aai_dir, 'shared_genes')
        files = os.listdir(shared_genes_dir)

        homologs = defaultdict(set)
        for f in files:
            with open(os.path.join(shared_genes_dir, f)) as fin:
                fin.readline()

                for line in fin:
                    line_split = line.split('\t')

                    gene_idA = line_split[0]
                    gene_idB = line_split[1]

                    homologs[gene_idA].add(gene_idB)
                    homologs[gene_idB].add(gene_idA)

        genes_to_remove = set()
        genes_to_keep = set()
        sorted_keys = sorted(homologs, key=lambda k: len(homologs[k]), reverse=True)
        for gene_id in sorted_keys:
            gene_set = homologs[gene_id]

            if len(gene_set.intersection(genes_to_keep)) > 0:
                genes_to_remove.update(gene_set - genes_to_keep)
                genes_to_remove.add(gene_id)
            else:
                genes_to_keep.add(gene_id)
                genes_to_remove.update(gene_set - genes_to_keep)

        # The CompareM call to rblast creates fasta files where gene ids are modified to
        # also contain genome ids. This is just a hack so to point to the directory with
        # these amended fasta files.
        os.system('ln -s %s %s' % (os.path.join(rblast_dir, 'genes'), ammended_gene_dir))

        return genes_to_remove

    def run(self,
                taxonomy_file, type_strains_file,
                genome_dir, max_taxa, rank,
                per_identity, per_aln_len,
                genomes_to_process, keep_all_genes,
                create_diamond_db, create_blast_db,
                cpus, output_dir):
        """ Create dereplicate set of genes.

        Taxonomy file should have the following format:
            <genome_id>\t<taxonomy_str>

            where taxonomy_str is in GreenGenes format:
                d__Bacteria;p__Proteobacteria;...;s__Escherichia coli

        Type strain file should have the following format:
            <genome_id>\t<genome name>

        Parameters
        ----------
        taxonomy_file : str
            File indicating taxonomy string for all genomes of interest
        type_strains_file : str
            File indicating type strains.
        genome_dir : str
            Directory with genomes in individual directories.
        max_taxa : int
            Maximum taxa to retain in a named group.
        rank : int
            Taxonomic rank to perform dereplication (0 = domain, ..., 6 = species).
        per_identity : float
            Percent identity for subsampling similar genes.
        per_aln_len : float
            Percent alignment length for subsampling similar genes.
        genomes_to_process : str
            File with list of genomes to retain instead of performing taxon subsampling.
        keep_all_genes : boolean
            Flag indicating that no gene subsampling should be performed.
        create_diamond_db : boolean
            Flag indicating if DIAMOND database should be created.
        create_blast_db : boolean
            Flag indicating if BLAST database should be created.
        cpus : int
            Number of cpus to use.
        output_dir : str
            Desired output directory for storing results.
        """

        make_sure_path_exists(output_dir)

        print 'Dereplicating at the rank of %s.' % self.rank_labels[rank]

        print ''
        print 'Reading taxonomy file.'
        taxonomy = self.read_taxonomy(taxonomy_file)
        print '  There are %d genomes with taxonomy strings.' % len(taxonomy)

        print ''
        print 'Reading type strain file.'
        type_strains = self.read_type_strain(type_strains_file)
        print '  There are %d type strains.' % len(type_strains)

        # get specific list of genomes to process
        genomes_to_retain = set()
        if genomes_to_process:
            print ''
            print 'Reading genomes to retain.'
            for line in open(genomes_to_process):
                line_split = line.split()
                genomes_to_retain.add(line_split[0])
            print '  Retaining %d genomes.' % len(genomes_to_retain)

        # identify unique genes in each named group
        fout = open(os.path.join(output_dir, 'genomes_without_called_genes.tsv'), 'w')
        rank_genomes = defaultdict(list)
        genomes_with_missing_data = set()
        underclassified_genomes = 0
        for genome_id, t in taxonomy.iteritems():
            if genomes_to_process and genome_id not in genomes_to_retain:
                continue

            genome_file = os.path.join(genome_dir, genome_id, genome_id + '.genes.faa')
            if not os.path.exists(genome_file):
                genomes_with_missing_data.add(genome_id)
                fout.write(genome_id + '\t' + ';'.join(taxonomy[genome_id]) + '\n')
                continue

            taxa = t[rank]
            if taxa[3:] == '':
                underclassified_genomes += 1
                rank_genomes[self.underclassified].append(genome_id)
            else:
                rank_genomes[taxa].append(genome_id)
        fout.close()

        total_genomes_to_process = sum([len(genome_list) for genome_list in rank_genomes.values()])

        print ''
        print 'Under-classified genomes automatically placed into the database: %d' % underclassified_genomes
        print 'Genomes with missing sequence data: %d' % len(genomes_with_missing_data)
        print ''
        print 'Total named groups: %d' % len(rank_genomes)
        print 'Total genomes to process: %d' % total_genomes_to_process

        # process each named group
        print ''
        gene_file = os.path.join(output_dir, 'genome_db.%s.genes.faa' % str(datetime.date.today()))
        gene_out = open(gene_file, 'w')

        taxonomy_out = open(os.path.join(output_dir, 'taxonomy.%s.tsv' % str(datetime.date.today())), 'w')

        tmp_dir = tempfile.mkdtemp()
        total_genes_removed = 0
        total_genes_kept = 0
        total_genomes_kept = 0
        processed_genomes = 0
        for taxa, genome_list in rank_genomes.iteritems():
            processed_genomes += len(genome_list)

            print ''
            print '-------------------------------------------------------------------------------'
            print ' Processing %s | Finished %d of %d (%.2f%%) genomes.' % (taxa, processed_genomes, total_genomes_to_process, processed_genomes * 100.0 / total_genomes_to_process)
            print self.time_keeper.get_time_stamp()
            print '-------------------------------------------------------------------------------'

            # create directory with selected genomes
            taxon_dir = os.path.join(tmp_dir, 'taxon')
            os.mkdir(taxon_dir)

            reduced_genome_list = genome_list
            if not genomes_to_process and taxa != self.underclassified:  # perform taxon subsampling
                reduced_genome_list = self.select_taxa(genome_list, taxonomy, type_strains, max_taxa)
            total_genomes_kept += len(reduced_genome_list)

            gene_dir = os.path.join(taxon_dir, 'genes')
            os.mkdir(gene_dir)
            for genome_id in reduced_genome_list:
                taxonomy_out.write(genome_id + '\t' + ';'.join(taxonomy[genome_id]) + '\n')
                cur_genome_dir = os.path.join(genome_dir, genome_id)
                self.img_gene_id_to_scaffold_id(cur_genome_dir, genome_id, gene_dir)

            # filter genes based on amino acid identity
            genes_to_remove = []
            amended_gene_dir = os.path.join(taxon_dir, 'ammended_genes')
            if keep_all_genes or taxa == self.underclassified:
                # modify gene identifiers to include genome ids
                self.amend_gene_identifies(gene_dir, amended_gene_dir)
            else:
                # filter genes on AAI
                genes_to_remove = self.filter_aai(taxon_dir, gene_dir, amended_gene_dir, per_identity, per_aln_len, cpus)

            print ''
            print '  Writing unique genes from genomes in %s.' % taxa
            genes_kept = self.write_gene_file(gene_out, amended_gene_dir, reduced_genome_list, taxonomy, genes_to_remove)

            print '    Retain %d of %d taxa.' % (len(reduced_genome_list), len(genome_list))
            print '    Genes to keep: %d' % genes_kept
            print '    Genes removed: %d' % len(genes_to_remove)

            total_genes_kept += genes_kept
            total_genes_removed += len(genes_to_remove)

            shutil.rmtree(taxon_dir)

        taxonomy_out.close()
        gene_out.close()

        print ''
        print 'Retain %d of %d (%.1f%%) genomes' % (total_genomes_kept, total_genomes_to_process, total_genomes_kept * 100.0 / (total_genomes_to_process))
        print '  Total genes kept: %d' % total_genes_kept
        print '  Total genes removed: %d (%.1f%%)' % (total_genes_removed, total_genes_removed * 100.0 / (total_genes_kept + total_genes_removed))

        if create_diamond_db:
            print ''
            print 'Creating DIAMOND database.'
            os.system('diamond makedb -b 10 -p 32 -d %s --in %s' % (gene_file, gene_file))
            print ''

        if create_blast_db:
            print ''
            print 'Creating BLAST database.'
            os.system('makeblastdb -dbtype prot -in %s' % gene_file)
            print ''

        shutil.rmtree(tmp_dir)

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_taxonomy', help='taxonomy for all reference genomes')
    parser.add_argument('type_strains', help='file specifying type strains that should not be filtered')
    parser.add_argument('genome_dir', help='directory containing reference genomes in individual folders')
    parser.add_argument('output_dir', help='directory to store results')
    parser.add_argument('-m', '--max_taxa', type=int, default=50, help='maximum taxa to retain in a named group')
    parser.add_argument('-r', '--rank', type=int, default=5, help='rank to preform dereplication [0=domain, 6=species]')
    parser.add_argument('-p', '--per_identity', type=float, default=90.0, help="percent identity for subsampling similar genes")
    parser.add_argument('-a', '--per_aln_len', type=float, default=90.0, help="percent alignment length for subsampling similar genes")
    parser.add_argument('--genomes_to_process', default=None, help='list of genomes to retain instead of performing taxon subsampling')
    parser.add_argument('--keep_all_genes', action="store_true", default=False, help='restricts filtering to taxa')
    parser.add_argument('--create_diamond_db', action="store_true", default=False, help='create DIAMOND database')
    parser.add_argument('--create_blast_db', action="store_true", default=False, help='create BLAST database')
    parser.add_argument('-c', '--cpus', type=int, default=32, help='number of cpus to use')

    args = parser.parse_args()

    try:
        makeDatabase = MakeDatabase()
        makeDatabase.run(args.input_taxonomy,
                         args.type_strains,
                         args.genome_dir,
                         args.max_taxa,
                         args.rank,
                         args.per_identity,
                         args.per_aln_len,
                         args.genomes_to_process,
                         args.keep_all_genes,
                         args.create_diamond_db,
                         args.create_blast_db,
                         args.cpus,
                         args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
