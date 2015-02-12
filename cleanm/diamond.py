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

import os
import sys
import subprocess
import logging
from collections import namedtuple


class Diamond(object):
    """Wrapper for running diamond."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger()

        self._check_for_diamond()

        self.cpus = cpus

    def homology_criteria(self, evalue, per_identity, max_target_seqs=1):
        """Criteria to use in homology search.

        Parameters
        ----------
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast.
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        """

        self.evalue = evalue
        self.per_identity = per_identity
        self.max_target_seqs = max_target_seqs

    def _check_for_diamond(self):
        """Check to see if BLAST is on the system before we try to run it."""
        try:
            subprocess.call(['diamond', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  Make sure diamond is on your system path.")
            sys.exit()

    def read_blast_table(self, table):
        """Generator function to read hits from a blast output table.

        The table should be a TSV file in blast format 6. This is
        also the format used by Diamond.

        Parameters
        ----------
        table : str
            Name of table to read.

        Yields
        ------
        namedtuple
            Information about blast hit.
        """

        BlastHit = namedtuple('BlastHit', """query_id
                                                subject_id
                                                perc_identity
                                                aln_length
                                                mismatch_count
                                                gap_open_count
                                                query_start
                                                query_end
                                                subject_start
                                                subject_end
                                                evalue
                                                bitscore""")

        for line in open(table):
            line_split = line.split('\t')
            hit = BlastHit(query_id=line_split[0],
                            subject_id=line_split[1],
                            perc_identity=float(line_split[2]),
                            aln_length=int(line_split[3]),
                            mismatch_count=int(line_split[4]),
                            gap_open_count=int(line_split[5]),
                            query_start=int(line_split[6]),
                            query_end=int(line_split[7]),
                            subject_start=int(line_split[8]),
                            subject_end=int(line_split[9]),
                            evalue=float(line_split[10]),
                            bitscore=float(line_split[11]))

            yield hit

    def blastx(self, nt_file, db_file, output_file):
        """Apply diamond blastx to a set of nucleotide sequences.

        Parameters
        ----------
        nt_file : str
            Fasta file with nucleotide sequences.
        db_file : str
            Diamond database.
        output_file : str
            File to store hits identified by diamond.
        """

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        self.logger.info('  Running diamond blastx with %d processes (be patient!)' % self.cpus)
        os.system('diamond blastx --compress 0 -p %d -q %s -d %s -e %f --id %f -k %d -o %s' % (self.cpus,
                                                                            nt_file,
                                                                            db_file,
                                                                            self.evalue,
                                                                            self.per_identity,
                                                                            self.max_target_seqs,
                                                                            output_file))


