###############################################################################
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
import logging
from collections import defaultdict

import biolib.seq_io as seq_io
from biolib.common import (make_sure_path_exists,
                           concatenate_files)
from biolib.blast_parser import BlastParser
from biolib.external.diamond import Diamond

from numpy import mean


class Reference(object):
    """Compare scaffolds to a specified set of reference genomes."""

    def __init__(self, cpus, output_dir):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        output_dir : str
            Directory to store results.
        """
        self.logger = logging.getLogger()

        self.cpus = cpus
        self.output_dir = output_dir

    def _top_hits_to_reference(self, hits_ref_genomes, hits_comp_ref_genomes):

        # get top hit to reference genomes of interest
        hits_to_ref = {}
        blast_parser = BlastParser()
        for hit in blast_parser.read_hit(hits_ref_genomes):
            if hit.query_id not in hits_to_ref:
                hits_to_ref[hit.query_id] = hit

        # see if there is a better hit to a competing genome
        for hit in blast_parser.read_hit(hits_comp_ref_genomes):
            if hit.query_id in hits_to_ref:
                if hit.bitscore > hits_to_ref[hit.query_id].bitscore:
                    del hits_to_ref[hit.query_id]

        return hits_to_ref

    def run(self, scaffold_gene_file, ref_genome_gene_files, db_file, evalue, per_identity, scaffold_id_bin_id):
        """Create taxonomic profiles for a set of genomes.

        Parameters
        ----------
        scaffold_gene_file : str
            Fasta file of genes residing of scaffolds
        ref_genome_gene_files : list of str
            Fasta files of called genes on reference genomes of interest.
        db_file : str
            Database of competing reference genes.
        evalue : float
            E-value threshold used by blast.
        per_identity: float
            Percent identity threshold used by blast.
        scaffold_id_bin_id : d[scaffold_id] -> bin_id
            Specifies id of binned scaffolds.
        """

        self.logger.info('')
        self.logger.info('  Creating diamond database for reference genomes.')
        ref_gene_file = os.path.join(self.output_dir, 'ref_genes.faa')
        concatenate_files(ref_genome_gene_files, ref_gene_file)

        diamond = Diamond(self.cpus)
        ref_diamond_db = os.path.join(self.output_dir, 'ref_genes')
        diamond.make_database(ref_gene_file, ref_diamond_db)

        self.logger.info('')
        self.logger.info('  Identifying homologs within reference genomes of interest.')
        self.diamond_dir = os.path.join(self.output_dir, 'diamond')
        make_sure_path_exists(self.diamond_dir)
        hits_ref_genomes_daa = os.path.join(self.diamond_dir, 'ref_hits')
        diamond.blastp(scaffold_gene_file, ref_diamond_db, evalue, per_identity, 1, hits_ref_genomes_daa)

        hits_ref_genomes = os.path.join(self.diamond_dir, 'ref_hits.tsv')
        diamond.view(hits_ref_genomes_daa + '.daa', hits_ref_genomes)

        self.logger.info('')
        self.logger.info('  Identifying homologs within competing reference genomes.')
        hits_comp_ref_genomes_daa = os.path.join(self.diamond_dir, 'competing_ref_hits')
        diamond.blastp(scaffold_gene_file, db_file, evalue, per_identity, 1, hits_comp_ref_genomes_daa)

        hits_comp_ref_genomes = os.path.join(self.diamond_dir, 'competing_ref_hits.tsv')
        diamond.view(hits_comp_ref_genomes_daa + '.daa', hits_comp_ref_genomes)

        # get list of genes with a top hit to the reference genomes of interest
        hits_to_ref = self._top_hits_to_reference(hits_ref_genomes, hits_comp_ref_genomes)

        # get number of genes on each scaffold
        num_genes_on_scaffold = defaultdict(int)
        for seq_id, _seq in seq_io.read_seq(scaffold_gene_file):
            scaffold_id = seq_id[0:seq_id.rfind('_')]
            num_genes_on_scaffold[scaffold_id] += 1

        # get hits to each scaffold
        hits_to_scaffold = defaultdict(list)
        for query_id, hit in hits_to_ref.iteritems():
            scaffold_id = query_id[0:query_id.rfind('_')]
            hits_to_scaffold[scaffold_id].append(hit)

        # report summary stats for each scaffold
        fout = open(os.path.join(self.output_dir, 'references.tsv'), 'w')
        fout.write('Scaffold Id\tBin Id\t# genes\t# hits\t% genes\tAvg. alignment length (bp)\tAvg. evalue\tAvg. bitscore\n')
        for scaffold_id, hits in hits_to_scaffold.iteritems():
            aln_len = []
            evalue = []
            bitscore = []
            for hit in hits:
                aln_len.append(hit.aln_length)
                evalue.append(hit.evalue)
                bitscore.append(hit.bitscore)
            fout.write('%s\t%s\t%d\t%d\t%.1f%%\t%d\t%.1g\t%.1f\n' % (scaffold_id,
                                                                    scaffold_id_bin_id.get(scaffold_id, 'unbinned'),
                                                                    num_genes_on_scaffold[scaffold_id],
                                                                    len(hits),
                                                                    len(hits) * 100.0 / num_genes_on_scaffold[scaffold_id],
                                                                    mean(aln_len),
                                                                    mean(evalue),
                                                                    mean(bitscore)))

        fout.close()
