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
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import remove_extension


class BinComparer(object):
    """Identify differences between two sets of genomes.

    Sequence within bins are identified by name. It is assumed
    genomes were constructed over a common set of contigs/scaffolds.
    """

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger('timestamp')
        self.reporter = logging.getLogger('no_timestamp')

    def _genome_seqs(self, genome_files):
        """Get unique id of sequences in each genome.

        Parameters
        ----------
        genome_files : iterable
            Genome files in fasta format.

        Returns
        -------
        dict: d[genome_id] -> set(seq_id1, ..., seq_idN)
            Ids of sequences in each genome.
        """

        genome_seqs = defaultdict(set)
        for genome_file in genome_files:
            genome_id = remove_extension(genome_file)
            for seq_id, _seq in seq_io.read_seq(genome_file):
                genome_seqs[genome_id].add(seq_id)

        return genome_seqs

    def _genome_stats(self, genomes, seq_lens):
        """Get basic statistics about genomes.

        Parameters
        ----------
        genomes : d[genome_id] ->  set(seq_id1, ..., seq_idN)
            Ids of sequences in each genome.
        seq_lens : iterable
            Length of sequences.

        Returns
        -------
        dict: d[genome_id] -> set(seq_id1, ..., seq_idN)
            Ids of sequences in each genome.
        """

        total_uniq_binned_seqs = 0
        total_uniq_binned_bases = 0

        genome_stats = {}
        processed_seqs = set()
        repeats = set()
        for binId, seqs in genomes.iteritems():
            num_binned_bases = 0
            for seq_id in seqs:
                num_binned_bases += seq_lens[seq_id]
                if seq_id not in processed_seqs:
                    processed_seqs.add(seq_id)
                    total_uniq_binned_bases += seq_lens[seq_id]
                    total_uniq_binned_seqs += 1
                else:
                    repeats.add(seq_id)

            genome_stats[binId] = [len(seqs), num_binned_bases]

        return genome_stats, total_uniq_binned_seqs, total_uniq_binned_bases, len(repeats)

    def run(self, genome_files1, genome_files2, seq_file, output_file):
        """Get basic statistics about genomes.

        Parameters
        ----------
        genome_files1 : iterable
            First set of henome files in fasta format.
        genome_files2 : iterable
            Second set of genome files in fasta format.
        seq_file : str
            Scaffolds/contigs binned to create genomes.
        output_file : str
            Desire file to write results.
        """

        # determine total number of sequences
        self.logger.info('Reading sequences.')

        seq_lens = {}
        total_bases = 0
        num_seqs_over_length = defaultdict(int)
        total_bases_over_length = defaultdict(int)
        lengths_to_check = [1000, 5000, 10000, 20000, 50000]
        for seq_id, seq in seq_io.read_seq(seq_file):
            seq_len = len(seq)
            seq_lens[seq_id] = seq_len
            total_bases += seq_len

            for length in lengths_to_check:
                if seq_len >= length:
                    num_seqs_over_length[length] += 1
                    total_bases_over_length[length] += seq_len

        # determine sequences in each bin
        genome_seqs1 = self._genome_seqs(genome_files1)
        genome_seqs2 = self._genome_seqs(genome_files2)

        # determine bin stats
        genome_stats1, total_uniq_binned_seqs1, total_uniq_binned_bases1, num_repeats1 = self._genome_stats(genome_seqs1, seq_lens)
        genome_stats2, total_uniq_binned_seqs2, total_uniq_binned_bases2, num_repeats2 = self._genome_stats(genome_seqs2, seq_lens)

        # sort bins by size
        genome_stats1 = sorted(genome_stats1.iteritems(), key=lambda x: x[1][1], reverse=True)
        genome_stats2 = sorted(genome_stats2.iteritems(), key=lambda x: x[1][1], reverse=True)

        # report summary results
        self.reporter.info('Total seqs = %d (%.2f Mbp)' % (len(seq_lens), float(total_bases) / 1e6))
        for length in lengths_to_check:
            self.reporter.info('  # seqs > %d kbp = %d (%.2f Mbp)' % (int(length / 1000),
                                                                        num_seqs_over_length[length],
                                                                        float(total_bases_over_length[length]) / 1e6))

        self.reporter.info('')
        self.reporter.info('Binned seqs statistics:')
        self.reporter.info('  1) # genomes: %s, # binned seqs: %d (%.2f%%), # binned bases: %.2f Mbp (%.2f%%), # seqs in multiple bins: %d'
                                % (len(genome_seqs1),
                                   total_uniq_binned_seqs2,
                                   float(total_uniq_binned_seqs1) * 100 / len(seq_lens),
                                   float(total_uniq_binned_bases1) / 1e6,
                                   float(total_uniq_binned_bases1) * 100 / total_bases,
                                   num_repeats1))
        self.reporter.info('  2) # genomes: %s, # binned seqs: %d (%.2f%%), # binned bases: %.2f Mbp (%.2f%%), # seqs in multiple bins: %d'
                                % (len(genome_seqs2),
                                   total_uniq_binned_seqs2,
                                   float(total_uniq_binned_seqs2) * 100 / len(seq_lens),
                                   float(total_uniq_binned_bases2) / 1e6,
                                   float(total_uniq_binned_bases2) * 100 / total_bases,
                                   num_repeats2))

        # output report
        fout = open(output_file, 'w')
        for data in genome_stats2:
            fout.write('\t' + data[0])
        fout.write('\tunbinned\t# seqs\t# bases (Mbp)\tBest match\t% bases in common\t% seqs in common\n')

        max_bp_common2 = defaultdict(int)
        max_seqs_common2 = defaultdict(int)
        best_matching_genome2 = {}
        binned_seqs2 = defaultdict(set)
        for data1 in genome_stats1:
            bin_id1 = data1[0]
            fout.write(bin_id1)

            seqs1 = genome_seqs1[bin_id1]

            max_bp_common = 0
            max_seqs_common = 0
            best_matching_genome = 'n/a'
            binned_seqs = set()
            for data2 in genome_stats2:
                bin_id2 = data2[0]
                seqs2 = genome_seqs2[bin_id2]

                seqs_common = seqs1.intersection(seqs2)
                binned_seqs.update(seqs_common)
                num_seqs_common = len(seqs_common)
                fout.write('\t' + str(num_seqs_common))

                bases_common = 0
                for seqId in seqs_common:
                    bases_common += seq_lens[seqId]

                if bases_common > max_bp_common:
                    max_bp_common = bases_common
                    max_seqs_common = num_seqs_common
                    best_matching_genome = bin_id2

                if bases_common > max_bp_common2[bin_id2]:
                    max_bp_common2[bin_id2] = bases_common
                    max_seqs_common2[bin_id2] = num_seqs_common
                    best_matching_genome2[bin_id2] = bin_id1

                binned_seqs2[bin_id2].update(seqs_common)
            fout.write('\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\n' % (len(seqs1) - len(binned_seqs),
                                                             data1[1][0],
                                                             float(data1[1][1]) / 1e6,
                                                             best_matching_genome,
                                                             float(max_bp_common) * 100 / data1[1][1],
                                                             float(max_seqs_common) * 100 / data1[1][0],
                                                             ))

        fout.write('unbinned')
        for data in genome_stats2:
            genome_id = data[0]
            fout.write('\t%d' % (len(genome_seqs2[genome_id]) - len(binned_seqs2[genome_id])))
        fout.write('\n')

        fout.write('# seqs')
        for data in genome_stats2:
            fout.write('\t%d' % data[1][0])
        fout.write('\n')

        fout.write('# bases (Mbp)')
        for data in genome_stats2:
            fout.write('\t%.2f' % (float(data[1][1]) / 1e6))
        fout.write('\n')

        fout.write('Best match')
        for data in genome_stats2:
            binId = data[0]
            fout.write('\t%s' % best_matching_genome2.get(binId, 'n/a'))
        fout.write('\n')

        fout.write('% bases in common')
        for data in genome_stats2:
            binId = data[0]
            fout.write('\t%.2f' % (float(max_bp_common2[binId]) * 100 / data[1][1]))
        fout.write('\n')

        fout.write('% seqs in common')
        for data in genome_stats2:
            binId = data[0]
            fout.write('\t%.2f' % (float(max_seqs_common2[binId]) * 100 / data[1][0]))
        fout.write('\n')

        fout.close()
