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

import sys
import logging
from collections import namedtuple, defaultdict

from refinem.coverage import Coverage
from refinem.tetranucleotide import Tetranucleotide
from refinem.errors import ParsingError

from biolib.common import remove_extension
import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk


"""
To Do:
 1. Should split run() method so it produces a dictionary of named tuples
        with scaffold statistics. A seperate function should be used to
        write this to file. This would mirror the interface of GenomeStats().
"""


class ScaffoldStats(object):
    """Statistics for scaffolds.

    This class holds general statistics for individual scaffolds:
     - bin assignment
     - gc
     - len
     - coverage
     - tetranucleotide signature
    """

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        self.unbinned = 'unbinned'

        self.ScaffoldStats = namedtuple('ScaffoldStats', """genome_id
                                                            gc
                                                            length
                                                            coverage
                                                            signature""")

    def run(self, scaffold_file, genome_files, tetra_file, coverage_file, output_file):
        """Calculate statistics for scaffolds.

        Parameters
        ----------
        scaffold_file : str
            Fasta file containing scaffolds.
        genome_files : list of str
            Fasta files with binned scaffolds.
        tetra_file : str
            Tetranucleotide signatures for scaffolds.
        coverage_file : str
            Coverage profiles for scaffolds
        output_file : str
            Output file for scaffolds statistics.
        """

        tetra = Tetranucleotide(self.cpus)
        signatures = tetra.read(tetra_file)

        cov_profiles = None
        if coverage_file:
            coverage = Coverage(self.cpus)
            cov_profiles, _ = coverage.read(coverage_file)

        # determine bin assignment for each scaffold
        self.logger.info('Determining scaffold statistics.')

        scaffold_id_genome_id = {}
        for gf in genome_files:
            genome_id = remove_extension(gf)
            for scaffold_id, _seq in seq_io.read_seq(gf):
                scaffold_id_genome_id[scaffold_id] = genome_id

        # write out scaffold statistics
        fout = open(output_file, 'w')
        fout.write('Scaffold id\tGenome Id\tGC\tLength (bp)')

        if cov_profiles:
            bam_ids = sorted(cov_profiles[cov_profiles.keys()[0]].keys())
            for bam_id in bam_ids:
                fout.write('\t' + bam_id)

        for kmer in tetra.canonical_order():
            fout.write('\t' + kmer)
        fout.write('\n')

        for scaffold_id, seq in seq_io.read_seq(scaffold_file):
            fout.write(scaffold_id)
            fout.write('\t' + scaffold_id_genome_id.get(scaffold_id, self.unbinned))
            fout.write('\t%.2f' % (seq_tk.gc(seq) * 100.0))
            fout.write('\t%d' % len(seq))

            if cov_profiles:
                for bam_id in bam_ids:
                    fout.write('\t%.2f' % cov_profiles[scaffold_id][bam_id])

            fout.write('\t' + '\t'.join(map(str, signatures[scaffold_id])))
            fout.write('\n')

        fout.close()

    def read(self, stats_file):
        """Read statistics for scaffolds.

        Parameters
        ----------
        stats_file : str
            File with statistics for individual scaffolds.
        """

        try:
            sig = {}
            self.genome_ids = set()
            with open(stats_file) as f:
                header = f.readline().split('\t')

                if 'AAAA' not in header:
                    raise ParsingError("[Error] Statistics file is missing tetranucleotide signature data: %s" % stats_file)

                tetra_index = header.index('AAAA')
                self.signature_headers = [x.strip() for x in header[tetra_index:]]
                self.coverage_headers = [x.strip() for x in header[4:tetra_index]]

                self.scaffolds_in_genome = defaultdict(set)
                self.stats = {}
                for line in f:
                    line_split = line.split('\t')
                    scaffold_id = line_split[0]
                    genome_id = line_split[1]
                    gc = float(line_split[2])
                    scaffold_len = int(line_split[3])

                    coverage = []
                    for cov in line_split[4:tetra_index]:
                        coverage.append(float(cov))

                    signature = []
                    for freq in line_split[tetra_index:]:
                        signature.append(float(freq))

                    self.stats[scaffold_id] = self.ScaffoldStats(genome_id, gc, scaffold_len, coverage, signature)

                    if genome_id != self.unbinned:
                        self.scaffolds_in_genome[genome_id].add(scaffold_id)

            return sig
        except IOError:
            print '[Error] Failed to open scaffold statistics file: %s' % stats_file
            sys.exit()
        except ParsingError:
            sys.exit()

    def num_scaffolds(self):
        """Number of scaffolds.

        Returns
        -------
        int
            Number of scaffolds.
        """

        return len(self.stats)

    def num_genomes(self):
        """Number of genomes.

        Returns
        -------
        int
            Number of genomes.
        """

        return len(self.scaffolds_in_genome)

    def coverage_profile_length(self):
        """Length of coverage profile.

        Returns
        -------
        int
            Length of coverage profile.
        """

        return len(self.coverage_headers)

    def signature_length(self):
        """Length of tetranucleotide signature.

        Returns
        -------
        int
            Length of tetranucleotide signature.
        """

        return len(self.signature_headers)

    def get(self, scaffold_id):
        """Statistics of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        namedtuple -> genome_id, gc, scaffold_len, coverage, signature
            Statistics for scaffold.
        """

        return self.stats[scaffold_id]

    def genome_id(self, scaffold_id):
        """Genome assignment of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        str
            Genome assignment of scaffold.
        """

        return self.stats[scaffold_id].genome_id

    def gc(self, scaffold_id):
        """GC of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        float
            GC of scaffold.
        """

        return self.stats[scaffold_id].gc

    def scaffold_length(self, scaffold_id):
        """Length of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        int
            Length of scaffold.
        """

        return self.stats[scaffold_id].length

    def coverage(self, scaffold_id):
        """Coverage profile of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        list
            Coverage profile of scaffold.
        """

        return self.stats[scaffold_id].coverage

    def signature(self, scaffold_id):
        """Tetranucleotide signature of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        list
           Tetranucleotide signature of scaffold.
        """

        return self.stats[scaffold_id].signature

    def print_coverage_header(self):
        """Print header line for coverage profile."""

        return '\t'.join(self.coverage_headers)

    def print_signature_header(self):
        """Print header line for tetranucleotide signature."""

        return '\t'.join(self.signature_headers)

    def print_stats(self, scaffold_id):
        """Produce string indicating binning, GC, and length of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        str
            String indicating genome id, scaffold length, and scaffold GC
        """

        stats = self.stats[scaffold_id]
        return '%s\t%d\t%.2f' % (stats.genome_id, stats.length, stats.gc)

    def print_coverage(self, scaffold_id):
        """Produce string indicating coverage profile of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        str
            Coverage profile for scaffold of interest.
        """

        cov_strs = []
        for cov in self.stats[scaffold_id].coverage:
            cov_strs.append('%.2f' % cov)

        return '\t'.join(cov_strs)

    def print_signature(self, scaffold_id):
        """Produce string indicating tetranucleotide signature of scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold of interest.

        Returns
        -------
        str
            Tetranucleotide signature for scaffold of interest.
        """

        tetra_strs = []
        for tetra in self.stats[scaffold_id].signature:
            tetra_strs.append('%.2f' % tetra)

        return '\t'.join(tetra_strs)
