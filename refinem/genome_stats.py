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

from numpy import (array as np_array, mean as np_mean, median as np_median)
import weightedstats as ws

from biolib.common import alphanumeric_sort
from biolib.genomic_signature import GenomicSignature


class GenomeStats():
    """Statistics for genomes.

    This class calculates statistics for genomes comprised
    of one or more scaffolds. Mean and median statistics are 
    weighted by scaffold length. The following weighted statistics 
    are calculated:
     - mean and median scaffold length
     - mean and median GC
     - mean and median coverage
     - mean and median tetranucleotide signature
     - mean and median tetranucleotide distance (TD) from mean/median of genome
    """

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.GenomeStats = namedtuple('GenomeStats', """genome_size
                                                        mean_scaffold_length
                                                        median_scaffold_length
                                                        mean_gc
                                                        median_gc
                                                        mean_coverage
                                                        median_coverage
                                                        mean_signature
                                                        mean_td
                                                        median_td
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

        genome_size = defaultdict(int)
        scaffold_length = defaultdict(list)
        gc = defaultdict(list)
        coverage = defaultdict(list)
        signature = defaultdict(list)
        for _scaffold_id, stats in scaffold_stats.stats.iteritems():
            if stats.genome_id == scaffold_stats.unbinned:
                continue
                
            genome_size[stats.genome_id] += stats.length
            scaffold_length[stats.genome_id].append(stats.length)
            gc[stats.genome_id].append(stats.gc)
            coverage[stats.genome_id].append(stats.coverage)
            signature[stats.genome_id].append(stats.signature)

        # record statistics for each genome
        genomic_signature = GenomicSignature(0)

        self.genome_stats = {}
        for genome_id in genome_size:
            # calculate weighted mean and median statistics
            weights = np_array(scaffold_length[genome_id])
            
            len_array = np_array(scaffold_length[genome_id])
            mean_len = ws.numpy_weighted_mean(len_array, weights)
            median_len = ws.numpy_weighted_median(len_array, weights)
            
            gc_array = np_array(gc[genome_id])
            mean_gc = ws.numpy_weighted_mean(gc_array, weights)
            median_gc = ws.numpy_weighted_median(gc_array, weights)
            
            cov_array = np_array(coverage[genome_id]).T
            mean_cov = ws.numpy_weighted_mean(cov_array, weights)
            median_cov = []
            for i in xrange(cov_array.shape[0]):
                median_cov.append(ws.numpy_weighted_median(cov_array[i,:], weights))
            
            signature_array = np_array(signature[genome_id]).T
            mean_signature = ws.numpy_weighted_mean(signature_array, weights)

            # calculate mean and median tetranucleotide distance
            td = []
            for scaffold_id in scaffold_stats.scaffolds_in_genome[genome_id]:
                stats = scaffold_stats.stats[scaffold_id]
                td.append(genomic_signature.manhattan(stats.signature, mean_signature))
  
            self.genome_stats[genome_id] = self.GenomeStats(genome_size[genome_id],
                                                            mean_len, median_len,
                                                            mean_gc, median_gc,
                                                            mean_cov, median_cov,
                                                            mean_signature,
                                                            np_mean(td), np_median(td))

        return self.genome_stats

    def write(self, output_file):
        """Write genome statistics to file.

        Parameters
        ----------
        output_file : str
            Name of output file.
        """

        fout = open(output_file, 'w')
        fout.write('Genome id\tGenome size (bp)')
        fout.write('\tMean GC\tMedian GC')
        fout.write('\tMean scaffold length (bp)\tMean scaffold length (bp)')
        fout.write('\tMean: ' + '\tMean: '.join(self.coverage_headers))
        fout.write('\tMedian: ' + '\tMedian: '.join(self.coverage_headers))
        fout.write('\t' + '\t'.join(self.signature_headers))
        fout.write('\n')

        for genome_id in alphanumeric_sort(self.genome_stats.keys()):
            stats = self.genome_stats[genome_id]

            fout.write(genome_id)
            fout.write('\t%d' % stats.genome_size)
            fout.write('\t%.2f' % stats.mean_gc)
            fout.write('\t%.2f' % stats.median_gc)
            fout.write('\t%.2f' % stats.mean_scaffold_length)
            fout.write('\t%.2f' % stats.median_scaffold_length)

            for cov in stats.mean_coverage:
                fout.write('\t%.2f' % cov)
                
            for cov in stats.median_coverage:
                fout.write('\t%.2f' % cov)

            for freq in stats.mean_signature:
                fout.write('\t%f' % freq)

            fout.write('\n')

        fout.close()
