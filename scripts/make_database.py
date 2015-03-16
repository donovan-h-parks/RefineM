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
import random
import tempfile
import shutil
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import make_sure_path_exists
from biolib.taxonomy import Taxonomy


class MakeDatabase(object):
    """Make a dereplicated database of genes.

    Dereplication is done between genes within a named taxonomic
    group (e.g., genomes in the same genus) and is based on the
    average amino acid (AAI) between genes. Groups with large
    numbers of taxa can take an excessive amount of time to
    process so can be subsampled to a specific number of taxa.
    Subsampling is done in a manor which aims to retain
    phylogenetic diversity and thus helps ensures a good
    distribution of genes within the group.
    """

    def __init__(self):
        """Initialize."""

        self.rank_prefixes = Taxonomy().rank_prefixes
        self.rank_index = Taxonomy().rank_index

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

    def select_taxa(self, genome_list, taxonomy, max_taxa):
        """Select subset of genomes with a good distribution across named groups.

        Groups genomes into named species groups and subsamples evenly across
        these groups. Some genomes may not have a species identifier. All such
        genome are simply treated as being from the same named species group
        (i.e., the 'unnamed' group).

        Parameters
        ----------
        genome_list : iterable
            Genome to subsample.
        taxonomy : d[genome_id] -> [domain, ..., species]
            Taxonomy of each genome.
        max_taxa : int
            Number of genomes to retain.

        Returns
        -------
        iterable
            Subsampled list of genomes.
        """

        if len(genome_list) <= max_taxa:
            return genome_list

        # group genomes into named species
        species = defaultdict(set)
        for genome_id in genome_list:
            taxa = taxonomy[genome_id][self.rank_index['s__']]
            species[taxa].add(genome_id)

        # sample genomes from each named species
        reduced_genome_list = []

        while len(reduced_genome_list) != max_taxa:
            genomes_to_select = max_taxa - len(reduced_genome_list)
            genomes_per_group = max(genomes_to_select / len(species), 1)
            for taxa, genome_ids in species.iteritems():
                selected_genomes = random.sample(genome_ids, min(len(genome_ids), genomes_per_group))
                species[taxa] = genome_ids.difference(selected_genomes)

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
        taxonomy : d[genome_id] -> [domain, ..., species]
            Taxonomy of each genome.
        genes_to_ignore : set
            Genes which should not be written to file.
        """

        for genome_id in genome_list:
            genome_gene_file = os.path.join(gene_dir, genome_id + '.genes.faa')
            if not os.path.exists(genome_gene_file):
                continue

            taxonomy_str = ';'.join(taxonomy[genome_id])

            if os.stat(genome_gene_file).st_size == 0:
                continue

            for gene_id, seq in seq_io.read_fasta_seq(genome_gene_file):
                if gene_id in genes_to_ignore:
                    continue

                gene_out.write('>' + gene_id + ' ' + taxonomy_str + '\n')
                gene_out.write(seq + '\n')

    def run(self, input_taxonomy, genome_dir, max_taxa, rank, cpus, output_dir):
        """ Create dereplicate set of genes.

        Taxonomy file should have the following format:
            <genome_id>\t<taxonomy_str>

            where taxonomy_str is in GreenGenes format:
                d__Bacteria;p__Firmicutes;...

        Parameters
        ----------
        input_taxonomy : str
            Taxonomy file.
        genome_dir : str
            Directory with genomes in individual directories.
        max_taxa : int
            Maximum taxa to retain in a named group.
        rank : int
            Taxonomic rank to perform dereplication (0 = domain, ..., 6 = species).
        cpus : int
            Number of cpus to use.
        output_dir : str
            Desired output directory for storing results.
        """

        make_sure_path_exists(output_dir)

        print 'Reading taxonomy file.'
        taxonomy = self.read_taxonomy(input_taxonomy)

        gene_file = os.path.join(output_dir, 'genome_db.genes.faa')
        gene_out = open(gene_file, 'w')

        # identify unique genes in each species
        rank_genomes = defaultdict(list)
        genomes_with_missing_data = set()
        underclassified_genomes = []
        for genome_id, t in taxonomy.iteritems():
            taxa = t[rank]
            if taxa[3:] == '':
                underclassified_genomes.append(genome_id)
                continue

            genome_file = os.path.join(genome_dir, genome_id, genome_id + '.fna')
            if not os.path.exists(genome_file):
                genomes_with_missing_data.add(genome_id)
                continue

            rank_genomes[taxa].append(genome_id)

        print ''
        print 'Underclassified genomes being ignored: %d' % len(underclassified_genomes)

        print ''
        print 'Genomes with missing sequence data: %d' % len(genomes_with_missing_data)

        tmp_dir = tempfile.mkdtemp()
        total_genes_removed = 0
        processed_genomes = 0
        for taxa, genome_list in rank_genomes.iteritems():
            processed_genomes += len(genome_list)

            print ''
            print '*******************************************************************************'
            print ' Processing %s | Finished %d of %d (%.2f%%) genomes.' % (taxa, processed_genomes, len(taxonomy), processed_genomes * 100.0 / len(taxonomy))
            print '*******************************************************************************'

            species_dir = os.path.join(tmp_dir, taxa.replace(' ', '_'))
            os.mkdir(species_dir)

            if len(genome_list) > max_taxa:
                genome_list = self.select_taxa(genome_list, taxonomy, max_taxa)

            for genome_id in genome_list:
                genome_file = os.path.abspath(os.path.join(genome_dir, genome_id, genome_id + '.fna'))
                os.system('ln -s %s %s' % (genome_file, os.path.join(species_dir, genome_id + '.fna')))

            os.system('comparem aai_wf -p 90 -a 90 -e 1e-10 -c %d %s %s' % (cpus, species_dir, os.path.join(species_dir, 'aai')))

            # identify homologs to be filtered
            print ''
            print '  Identifying homologous genes to be filtered.'
            shared_genes_dir = os.path.join(species_dir, 'aai', 'shared_genes')
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

                if gene_id in genes_to_keep:
                    genes_to_remove.update(gene_set - genes_to_keep)
                elif len(gene_set.intersection(genes_to_keep)) > 0:
                    genes_to_remove.update(gene_set - genes_to_keep)
                    genes_to_remove.add(gene_id)
                else:
                    genes_to_keep.add(gene_id)
                    genes_to_remove.update(gene_set - genes_to_keep)

            print '    Number of taxa: %d' % len(genome_list)
            print '    Genes to keep: %d' % len(genes_to_keep)
            print '    Genes removed: %d' % len(genes_to_remove)

            total_genes_removed += len(genes_to_remove)

            print ''
            print '  Writing unique genes from genomes in %s.' % taxa
            gene_dir = os.path.join(species_dir, 'aai', 'genes')
            self.write_gene_file(gene_out, gene_dir, genome_list, taxonomy, genes_to_remove)

            shutil.rmtree(species_dir)

        print ''
        print 'Total genes removed: %d' % total_genes_removed

        print ''
        print 'Creating diamond database.'
        os.system('diamond makedb -b 10 -p 32 -d %s --in %s' % (gene_file, gene_file))
        print ''

        shutil.rmtree(tmp_dir)

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_taxonomy', help='taxonomy for all reference genomes')
    parser.add_argument('genome_dir', help='directory containing reference genomes in individual folders')
    parser.add_argument('output_dir', help='directory to store results')
    parser.add_argument('-m', '--max_taxa', type=int, default=50, help='maximum taxa to retain in a named group')
    parser.add_argument('-r', '--rank', type=int, default=5, help='rank to preform dereplication [0=domain, 6=species]')
    parser.add_argument('-c', '--cpus', type=int, default=32, help='number of cpus to use')

    args = parser.parse_args()

    try:
        makeDatabase = MakeDatabase()
        makeDatabase.run(args.input_taxonomy, args.genome_dir, args.max_taxa, args.rank, args.cpus, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
