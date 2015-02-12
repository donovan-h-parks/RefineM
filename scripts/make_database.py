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
__prog_desc__ = 'create database of reference genomes'

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
from collections import defaultdict


class MakeDatabase(object):
    def __init__(self):
        pass

    def read_taxonomy(self, input_taxonomy):
        taxonomy = {}
        for line in open(input_taxonomy):
            line_split = line.split('\t')
            taxonomy[line_split[0]] = [x.strip() for x in line_split[1].split(';')]

        return taxonomy

    def sample(self, taxonomy, max_species):
        # get all public genomes with a defined genus
        species = defaultdict(list)
        genomes_for_db = set()
        for genome_id, taxonomy in taxonomy.iteritems():
            if genome_id[0] == 'C':
                # ignore private genomes
                continue

            if taxonomy[5] == 'g__':
                # ignore genomes without a defined genus
                continue

            # tabulate genomes from the same named species
            if taxonomy[6] != 's__':
                species[taxonomy[6]].append(genome_id)
            else:
                genomes_for_db.add(genome_id)

        # subsample species with multiple representatives
        for genome_ids in species.values():
            if len(genome_ids) > max_species:
                genome_ids = random.sample(genome_ids, max_species)

            genomes_for_db.update(genome_ids)

        return genomes_for_db

    def write_taxonomy(self, genome_ids, taxonomy, taxonomy_file):
        fout = open(taxonomy_file, 'w')
        for genome_id in genome_ids:
            fout.write(genome_id + '\t' + ';'.join(taxonomy[genome_id]) + '\n')
        fout.close()

    def write_gene_file(self, genome_ids, genome_dir, gene_file):
        missing_genomes = set()
        with open(gene_file, 'w') as fout:
            for genome_id in genome_ids:
                genome_gene_file = os.path.join(genome_dir, genome_id + '.genes.faa')

                if not os.path.exists(genome_gene_file):
                    missing_genomes.add(genome_id)
                    continue

                gene_count = 0
                with open(genome_gene_file) as f:
                    for line in f:
                      if line[0] == '>':
                        fout.write('>' + genome_id + '~' + str(gene_count) + '\n')
                        gene_count += 1
                      else:
                        fout.write(line)

        return missing_genomes

    def run(self, input_taxonomy, genome_dir, max_species, output_dir):
        print 'Reading taxonomy file.'
        taxonomy = self.read_taxonomy(input_taxonomy)

        print 'Sampling %d genomes for each named species.' % max_species
        genomes_for_db = self.sample(taxonomy, max_species)
        print '  Identified %d reference genomes.' % len(genomes_for_db)

        print 'Writing genes from genomes.'
        gene_file = os.path.join(output_dir, 'genome_db.genes.faa')
        missing_genomes = self.write_gene_file(genomes_for_db, genome_dir, gene_file)
        genomes_for_db = genomes_for_db - missing_genomes
        print '  Removing %d reference genomes with no gene file.' % len(missing_genomes)

        print 'Writing taxonomy file.'
        taxonomy_file = os.path.join(output_dir, 'taxonomy.tsv')
        self.write_taxonomy(genomes_for_db, taxonomy, taxonomy_file)

        print 'Creating diamond database.'
        os.system('diamond makedb -b 10 -p 32 -d %s --in %s' % (gene_file, gene_file))

if __name__ == '__main__':
    print __prog_name__ + ' v' + __version__ + ': ' + __prog_desc__
    print '  by ' + __author__ + ' (' + __email__ + ')' + '\n'

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_taxonomy', help='taxonomy for all reference genomes')
    parser.add_argument('genome_dir', help='directory containing called genes for reference genomes')
    parser.add_argument('output_dir', help='directory to store results')

    parser.add_argument('-m', '--max_species', type=int, default=10, help='maximum genomes to sample for a named species')

    args = parser.parse_args()

    try:
        makeDatabase = MakeDatabase()
        makeDatabase.run(args.input_taxonomy, args.genome_dir, args.max_species, args.output_dir)
    except SystemExit:
        print "\nControlled exit resulting from an unrecoverable error or warning."
    except:
        print "\nUnexpected error:", sys.exc_info()[0]
        raise
