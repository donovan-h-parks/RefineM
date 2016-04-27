###############################################################################
#
# binStatistics.py - calculate statistics for each putative genome bin
#
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

import logging
from collections import namedtuple, defaultdict

from numpy import mean

from biolib.common import alphanumeric_sort
from biolib.genomic_signature import GenomicSignature


class GenomeStats():
    """Statistics for genomes.

    This class calculates statistics for genomes comprised
    of one or more scaffolds. Mean statistics are weighted by
    scaffold length. The following statistics are calculated:
     - bin assignment
     - mean GC
     - mean scaffold length
     - mean coverage
     - mean tetranucleotide signature
     - mean tetranucleotide distance (TD) from mean of genome
    """

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.GenomeStats = namedtuple('GenomeStats', """genome_size
                                                        mean_gc
                                                        mean_scaffold_length
                                                        mean_coverage
                                                        mean_signature
                                                        mean_td
                                                        """)

    def run(self, scaffold_stats):
        """Calculate statistics for genomes.

        Parameters
        ----------
        scaffold_stats : ScaffoldStats
            Statistics for individual scaffolds.
        """

        self.logger.info("Calculating statistics for %d genomes over %d scaffolds." % (scaffold_stats.num_genomes(),
                                                                                        scaffold_stats.num_scaffolds()))

        self.coverage_headers = scaffold_stats.coverage_headers
        self.signature_headers = scaffold_stats.signature_headers

        running_genomes_stats = defaultdict(lambda:
                                            [0,  # genome size
                                             0,  # mean_gc
                                             0,  # mean scaffold length
                                             [0] * scaffold_stats.coverage_profile_length(),  # coverage profile
                                             [0] * scaffold_stats.signature_length()])  # genomic signature

        for _scaffold_id, stats in scaffold_stats.stats.iteritems():
            if stats.genome_id == scaffold_stats.unbinned:
                continue

            genome_size, mean_gc, mean_length, mean_coverage, mean_signature = running_genomes_stats[stats.genome_id]

            genome_size += stats.length
            weight = float(stats.length) / genome_size

            mean_gc = stats.gc * weight + mean_gc * (1.0 - weight)
            mean_length = stats.length * weight + mean_length * (1.0 - weight)

            for i, cov in enumerate(stats.coverage):
                mean_coverage[i] = cov * weight + mean_coverage[i] * (1.0 - weight)

            for i, freq in enumerate(stats.signature):
                mean_signature[i] = freq * weight + mean_signature[i] * (1.0 - weight)

            running_genomes_stats[stats.genome_id] = genome_size, mean_gc, mean_length, mean_coverage, mean_signature

        # record statistics for each genome
        genomic_signature = GenomicSignature(0)

        self.genome_stats = {}
        for genome_id, gs in running_genomes_stats.iteritems():
            # calculate mean tetranucleotide distance
            td = []
            for scaffold_id in scaffold_stats.scaffolds_in_genome[genome_id]:
                stats = scaffold_stats.stats[scaffold_id]
                td.append(genomic_signature.manhattan(stats.signature, gs[4]))

            self.genome_stats[genome_id] = self.GenomeStats(gs[0], gs[1], gs[2], gs[3], gs[4], mean(td))

        return self.genome_stats

    def write(self, output_file):
        """Write genome statistics to file.

        Parameters
        ----------
        output_file : str
            Name of output file.
        """

        fout = open(output_file, 'w')
        fout.write('Genome id\tGenome size (bp)\tMean GC\tMean scaffold length (bp)')
        fout.write('\t' + '\t'.join(self.coverage_headers))
        fout.write('\t' + '\t'.join(self.signature_headers))
        fout.write('\n')

        for genome_id in alphanumeric_sort(self.genome_stats.keys()):
            stats = self.genome_stats[genome_id]

            fout.write(genome_id)
            fout.write('\t%d' % stats.genome_size)
            fout.write('\t%.2f' % stats.mean_gc)
            fout.write('\t%.2f' % stats.mean_scaffold_length)

            for cov in stats.mean_coverage:
                fout.write('\t%.2f' % cov)

            for freq in stats.mean_signature:
                fout.write('\t%f' % freq)

            fout.write('\n')

        fout.close()
