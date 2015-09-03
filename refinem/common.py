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

from biolib.common import remove_extension

import biolib.seq_io as seq_io


def concatenate_gene_files(gene_files, concatenated_gene_file):
    """Combine all gene files into a single file.

    Gene ids are modified to include genome ids in order to ensure
    all gene identifiers are unique across the set of genomes.

    Parameters
    ----------
    gene_files : list of str
        Fasta files of called genes to process.
    concatenated_gene_file : str
        Name of file to contain concatenated gene files.
    """

    fout = open(concatenated_gene_file, 'w')

    for gf in gene_files:
        genome_id = remove_extension(gf)

        for seq_id, seq in seq_io.read_seq(gf):
            fout.write('>' + seq_id + '~' + genome_id + '\n')
            if seq[-1] == '*':
                seq = seq[0:-1]
            fout.write(seq + '\n')

    fout.close()
