###############################################################################
#
# coverage.py - calculate coverage of all sequences
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

import os
import sys
import ast
import itertools
import logging

from scipy.stats import pearsonr
from numpy import (mean as np_mean)

import biolib.seq_io as seq_io
from biolib.common import find_nearest
from biolib.genomic_signature import GenomicSignature


class Outliers():
    """Identify scaffolds with divergent genomic characteristics."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

    def remove(self, genome_file, outlier_file, out_genome):
        """Remove sequences specified as outliers.

        Any scaffolds lists in the first column of
        the outlier file are removed from the specified
        genome.

        Parameters
        ----------
        genome_file : str
            Fasta file of binned scaffolds.
        outlier_file : str
            File specifying outlying scaffolds.
        out_genome : str
            Name of output genome.
        """

        genome_seqs = seq_io.read(genome_file)

        # remove scaffolds
        with open(outlier_file) as f:
            f.readline()

            for line in f:
                line_split = line.split('\t')
                scaffold_id = line_split[0]
                genome_seqs.pop(scaffold_id, None)

        # save modified bin
        seq_io.write_fasta(genome_seqs, out_genome)

    def identify(self, scaffold_stats, genome_stats,
                        gc_per, td_per,
                        cov_corr, cov_perc,
                        report_type, output_file):
        """Identify scaffolds with divergent genomic characteristics.

        Outliers are identified independently based on GC content,
        tetranucleotide signatures, coverage profile correlation, and
        mean percent deviation of coverage profile. The coverage correlation
        check is ignored if the coverage profile consists of a single value.

        Parameters
        ----------
        scaffold_stats : ScaffoldStats
            Statistics for individual scaffolds
        genome_stats : GenomeStats
            Staitsics for individual genomes.
        gc_per : int
            Percentile for identifying GC outliers
        td_per : int
            Percentile for identifying TD outliers
        cov_corr : int
            Correlation for identifying divergent coverage profiles
        cov_perc : int
            Mean percent deviation for identifying divergent coverage profiles
        report_type : str
            Report scaffolds that are outliers in 'all' or 'any' distribution
        output_file : str
            Name of output file.
        """

        # read reference distributions from file
        self.logger.info('  Reading reference distributions.')
        self.gc_dist = self._read_distribution('gc_dist')
        self.td_dist = self._read_distribution('td_dist')

        # identify outliers in each genome
        fout = open(output_file, 'w')
        fout.write('Scaffold id\tGenome id\tScaffold length (bp)\tOutlying distributions')
        fout.write('\tScaffold GC\tMean genome GC\tLower GC bound (%s%%)\tUpper GC bound (%s%%)' % (gc_per, gc_per))
        fout.write('\tScaffold TD\tMean genome TD\tUpper TD bound (%s%%)' % td_per)
        fout.write('\tMean scaffold coverage\tMean genome coverage\tCoverage correlation\tMean coverage deviation\n')

        genomic_signature = GenomicSignature(0)

        processed_genomes = 0
        for genome_id, scaffold_ids in scaffold_stats.scaffolds_in_genome.iteritems():
            processed_genomes += 1

            sys.stdout.write('    Finding outliers in %d of %d (%.1f%%) genomes.\r' % (processed_genomes,
                                                                                     scaffold_stats.num_genomes(),
                                                                                     processed_genomes * 100.0 / scaffold_stats.num_genomes()))
            sys.stdout.flush()

            # find keys into GC and TD distributions
            # gc -> [mean GC][scaffold length][percentile]
            # td -> [scaffold length][percentile]
            gs = genome_stats[genome_id]
            closest_gc = find_nearest(self.gc_dist.keys(), gs.mean_gc / 100.0)
            sample_seq_len = self.gc_dist[closest_gc].keys()[0]
            d = self.gc_dist[closest_gc][sample_seq_len]
            gc_lower_bound_key = find_nearest(d.keys(), (100 - gc_per) / 2.0)
            gc_upper_bound_key = find_nearest(d.keys(), (100 + gc_per) / 2.0)

            td_bound_key = find_nearest(self.td_dist[self.td_dist.keys()[0]].keys(), td_per)

            for scaffold_id in scaffold_ids:
                stats = scaffold_stats.stats[scaffold_id]

                # find GC and TD bounds
                closest_seq_len = find_nearest(self.gc_dist[closest_gc].keys(), stats.length)
                gc_lower_bound = self.gc_dist[closest_gc][closest_seq_len][gc_lower_bound_key]
                gc_upper_bound = self.gc_dist[closest_gc][closest_seq_len][gc_upper_bound_key]

                closest_seq_len = find_nearest(self.td_dist.keys(), stats.length)
                td_bound = self.td_dist[closest_seq_len][td_bound_key]

                # find changes from mean
                delta_gc = (stats.gc - gs.mean_gc) / 100.0
                delta_td = genomic_signature.manhattan(stats.signature, gs.mean_signature)

                # determine if scaffold is an outlier
                outlying_dists = []
                if delta_gc < gc_lower_bound or delta_gc > gc_upper_bound:
                    outlying_dists.append('GC')

                if delta_td > td_bound:
                    outlying_dists.append('TD')

                corr_r = 1.0
                if len(gs.mean_coverage) > 1:
                    corr_r, _corr_p = pearsonr(gs.mean_coverage, stats.coverage)
                    if  corr_r < cov_corr:
                        outlying_dists.append('COV_CORR')

                mean_cp = []
                for cov_genome, cov_scaffold in itertools.izip(gs.mean_coverage, stats.coverage):
                    mean_cp.append(abs(cov_genome - cov_scaffold) * 100 / cov_genome)

                mean_cp = np_mean(mean_cp)
                if mean_cp > cov_perc:
                    outlying_dists.append('COV_PERC')

                # report outliers
                if (report_type == 'any' and len(outlying_dists) >= 1) or (report_type == 'all' and len(outlying_dists) == 2):
                    fout.write('%s\t%s\t%s\t%s' % (scaffold_id, genome_id, stats.length, ','.join(outlying_dists)))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (stats.gc, gs.mean_gc, gs.mean_gc + gc_lower_bound * 100, gs.mean_gc + gc_upper_bound * 100))
                    fout.write('\t%.3f\t%.3f\t%.3f' % (delta_td, gs.mean_td, td_bound))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (np_mean(stats.coverage), np_mean(gs.mean_coverage), corr_r, mean_cp))
                    fout.write('\n')

        sys.stdout.write('\n')
        fout.close()

    def _read_distribution(self, prefix):
        """Read distribution file.

        Parameters
        ----------
        prefix : str
            Prefix of distibution to read (gc_dist or td_dist)

        Returns
        -------
        dict : d[percentile] -> critical value
            Critical values at integer percentiles.
        """

        dist_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'distributions', prefix + '.txt')

        with open(dist_file, 'r') as f:
            s = f.read()
            d = ast.literal_eval(s)

        return d
