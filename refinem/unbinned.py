###############################################################################
#
# unbinned.py - identify unbinned sequences
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

from biolib.common import check_file_exists

import biolib.seq_io as seq_io


class Unbinned():
    """Identified scaffolds not assigned to a putative genome."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger('timestamp')

    def run(self, genome_files, scaffold_file, min_seq_len):
        """Fragment genome sequences into fragments of a fixed size.

        Parameters
        ----------
        genome_files : list of str
            Fasta files of genomes to process.
        scaffold_file : str
            Scaffolds binned to generate putative genomes.
        min_seq_len : int
            Ignore scaffolds shorter than the specified length.

        Returns
        -------
        dict : d[seq_id] -> seq
            Dictionary of unbinned sequences.
        """

        check_file_exists(scaffold_file)

        # get list of sequences in bins
        self.logger.info('Reading binned scaffolds.')

        binned_seq_ids = set()
        total_binned_bases = 0
        for genome_file in genome_files:
            for seq_id, seq in seq_io.read_seq(genome_file):
                binned_seq_ids.add(seq_id)
                total_binned_bases += len(seq)

        self.logger.info('Read %d (%.2f Mbp) binned scaffolds.' % (len(binned_seq_ids), float(total_binned_bases) / 1e6))

        # write all unbinned sequences
        self.logger.info('Identifying unbinned scaffolds >= %d bp.' % min_seq_len)

        unbinned_bases = 0
        unbinned_seqs = {}
        for seq_id, seq in seq_io.read_seq(scaffold_file):
            if seq_id not in binned_seq_ids and len(seq) >= min_seq_len:
                unbinned_seqs[seq_id] = seq
                unbinned_bases += len(seq)

        self.logger.info('Identified %d (%.2f Mbp) unbinned scaffolds.' % (len(unbinned_seqs), float(unbinned_bases) / 1e6))

        self.logger.info('Percentage of unbinned scaffolds: %.2f%%' % (len(unbinned_seqs) * 100.0 / (len(unbinned_seqs) + len(binned_seq_ids))))
        self.logger.info('Percentage of unbinned bases: %.2f%%' % (unbinned_bases * 100.0 / (unbinned_bases + total_binned_bases)))

        return unbinned_seqs
