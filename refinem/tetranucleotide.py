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

from biolib.genomic_signature import GenomicSignature
from biolib.parallel import Parallel

import numpy as np

from refinem.errors import ParsingError


class Tetranucleotide(object):
    """Calculate tetranucleotide signature of sequences."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        self.signatures = GenomicSignature(4)

    def canonical_order(self):
        """Canonical order of tetranucleotides."""
        return self.signatures.canonical_order()

    def _producer(self, seq_info):
        """Calculate tetranucleotide signature of a sequence.

        Parameters
        ----------
        seq_id : str
            Unique id of sequence.
        seq : str
            Sequence in nuceltoide space.

        Returns
        -------
        str
            Unique id of sequence.
        list
            Count of each kmer in the canonical order.
        """

        seq_id, seq = seq_info

        sig = self.signatures.seq_signature(seq)

        total_kmers = sum(sig)
        for i in xrange(0, len(sig)):
            sig[i] = float(sig[i]) / total_kmers

        return (seq_id, sig)

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : list -> kmers in the canonical order
            Tetranucleotide signature in canconical order.
        consumer_data : d[seq_id] -> tetranucleotide signature
            Set of kmers observed across all genomes (kmer_set),
            along with the kmer usage of each genome (genome_kmer_usage).

        Returns
        -------
        consumer_data: dict
            The consumer data structure or None must be returned
        """

        if consumer_data == None:
            consumer_data = {}

        seq_id, sig = produced_data
        consumer_data[seq_id] = sig

        return consumer_data

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of sequences processed.
        total_items : int
            Total number of sequences to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        if self.logger.is_silent:
            return None
        else:
            return '  Finished processing %d of %d (%.2f%%) sequences.' % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self, seq_file):
        """Calculate tetranucleotide signatures of sequences.

        Parameters
        ----------
        seq_file : str
            Name of fasta/q file to read.

        Returns
        -------
        dict : d[seq_id] -> tetranucleotide signature in canonical order
            Count of each kmer.
        """

        self.logger.info('Calculating tetranucleotide signature for each sequence:')

        parallel = Parallel(self.cpus)
        seq_signatures = parallel.run_seqs_file(self._producer, self._consumer, seq_file, self._progress)

        return seq_signatures

    def read(self, signature_file):
        """Read tetranucleotide signatures.

        Parameters
        ----------
        signature_file : str
            Name of file to read.

        Returns
        -------
        dict : d[seq_id] -> tetranucleotide signature in canonical order
            Count of each kmer.
        """

        try:
            sig = {}
            with open(signature_file) as f:
                header = f.readline().split('\t')
                kmer_order = [x.strip().upper() for x in header[1:]]
                if len(kmer_order) != len(self.canonical_order()):
                    raise ParsingError("[Error] Tetranucleotide file must contain exactly %d tetranucleotide columns." % len(self.canonical_order()))

                canonical_order_index = np.argsort(kmer_order)
                canonical_order = [kmer_order[i] for i in canonical_order_index]

                if canonical_order != self.canonical_order():
                    raise ParsingError("[Error] Failed to process tetranucleotide signature file: " + signature_file)

                for line in f:
                    line_split = line.split('\t')
                    sig[line_split[0]] = [float(line_split[i + 1]) for i in canonical_order_index]

            return sig
        except IOError:
            print '[Error] Failed to open signature file: %s' % signature_file
            sys.exit()
        except ParsingError:
            sys.exit()

    def write(self, signatures, output_file):
        """Write tetranucleotide signatures.

        Parameters
        ----------
        signature_file : d[seq_id] -> tetranucleotide signature in canonical order
            Count of each kmer.
        output_file : str
            Name of output file.
        """

        fout = open(output_file, 'w')

        fout.write('Scaffold id')
        for kmer in self.canonical_order():
            fout.write('\t' + kmer)
        fout.write('\n')

        for seq_id, tetra_signature in signatures.iteritems():
            fout.write(seq_id + '\t')
            fout.write('\t'.join(map(str, tetra_signature)))
            fout.write('\n')

        fout.close()
