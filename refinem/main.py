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
import logging

from refinem.taxonomic_profile import TaxonomicProfile
from refinem.gene_profile import GeneProfile
from refinem.bin_comparer import BinComparer
from refinem.reference import Reference
from refinem.unbinned import Unbinned
from refinem.coverage import Coverage
from refinem.tetranucleotide import Tetranucleotide

import biolib.seq_io as seq_io
import biolib.genome_tk as genome_tk
from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           remove_extension,
                           check_file_exists,
                           query_yes_no)
from biolib.external.prodigal import Prodigal
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.execute import check_dependencies


"""
To do:
1. Need to consider when a command should take called genes instead of nucleotide sequences.
2. The reference command should only call genes on scaffolds that are not binned, otherwise it
   should use the called genes for the bin.
"""


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def _genome_files(self, genome_dir, genome_ext):
        """Identify genomes files.

        Parameters
        ----------
        genome_dir : str
            Directory containing genomes of interest.
        genome_ext : str
            Extension of genome files.

        Returns
        -------
        list
            Path to genome files.
        """

        check_dir_exists(genome_dir)

        genome_files = []
        for f in os.listdir(genome_dir):
            if f.endswith(genome_ext):
                genome_files.append(os.path.join(genome_dir, f))

        if not genome_files:
            self.logger.warning('  [Warning] No genomes found. Check the --genome_ext flag used to identify genomes.')
            sys.exit()

        return genome_files

    def _check_nuclotide_seqs(self, seq_files):
        """Check if files contain sequences in nucleotide space.

        Parameters
        ----------
        seq_files : iterable
            Sequence files to check.

        Returns
        -------
        boolean
            True if files can be treated as containing nucleotide sequences.
        """

        for seq_file in seq_files:
            if not seq_io.is_nucleotide(seq_file):
                print('Expected all files to contain sequenes in nucleotide space.')
                print('File %s appears like it may contain amino acids sequences.' % seq_file)

                yes_response = query_yes_no('Do all files contain only nucleotide sequences?', default='no')
                if not yes_response:
                    return False

        return True

    def taxa_profile(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - taxa_profile] Generating taxonomic profiles from short fragments.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)
        check_file_exists(options.taxonomy_file)
        check_file_exists(options.db_file)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        taxonomic_profile = TaxonomicProfile(options.cpus, options.output_dir)
        taxonomic_profile.run(genome_files,
                                 options.db_file,
                                 options.taxonomy_file,
                                 options.evalue,
                                 options.per_identity,
                                 options.window_size,
                                 options.step_size)

        self.logger.info('')
        self.logger.info('  Results written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def gene_profile(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - gene_profile] Generating taxonomic profiles from genes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)
        check_file_exists(options.taxonomy_file)
        check_file_exists(options.db_file)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        # call genes for each genome
        prodigal = Prodigal(options.cpus)
        gene_dir = os.path.join(options.output_dir, 'genes')
        make_sure_path_exists(gene_dir)
        prodigal.run(genome_files, False, None, False, gene_dir)

        # modify gene ids to include genome ids in order to ensure
        # all gene identifiers are unique across the set of genomes,
        # also removes the trailing asterisk used to identify the stop
        # codon
        self.logger.info('')
        self.logger.info('  Appending genome identifiers to all gene identifiers.')
        aa_gene_files = []
        for gf in genome_files:
            genome_id = remove_extension(gf)

            nt_file = os.path.join(gene_dir, genome_id + '.genes.fna')
            seqs = seq_io.read(nt_file)
            fout = open(nt_file, 'w')
            for seq_id, seq in seqs.iteritems():
                fout.write('>' + seq_id + '~' + genome_id + '\n')
                fout.write(seq + '\n')
            fout.close()

            aa_file = os.path.join(gene_dir, genome_id + '.genes.faa')
            seqs = seq_io.read(aa_file)
            fout = open(aa_file, 'w')
            for seq_id, seq in seqs.iteritems():
                fout.write('>' + seq_id + '~' + genome_id + '\n')
                if seq[-1] == '*':
                    seq = seq[0:-1]
                fout.write(seq + '\n')
            fout.close()

            aa_gene_files.append(aa_file)

        # build gene profile
        gene_profile = GeneProfile(options.cpus, options.output_dir)
        gene_profile.run(aa_gene_files,
                             options.db_file,
                             options.taxonomy_file,
                             options.evalue,
                             options.per_identity)

        self.logger.info('')
        self.logger.info('  Results written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def modify(self, options):
        """Modify command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - modify] Modifying scaffolds in genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(os.path.dirname(options.output_genome))

        if not (options.add or options.remove or options.outlier_file):
            self.logger.warning('  [Warning] No modification to bin requested.\n')
            sys.exit()

        if (options.add or options.remove) and options.outlier_file:
            self.logger.warning("  [Warning] The 'outlier_file' option cannot be specified with 'add' or 'remove'.\n")
            sys.exit()

        if options.add or options.remove:
            failed_to_add, failed_to_remove = genome_tk.modify(options.genome_file,
                                                               options.scaffold_file,
                                                               options.add,
                                                               options.remove,
                                                               options.output_genome)
        elif options.outlier_file:
            # **** TO DO!
            # binTools.removeOutliers(options.genome_file, options.outlier_file, options.output_file)
            pass

        if failed_to_add:
            self.logger.warning('  [Warning] Failed to add the following sequence(s):')
            for seq_id in failed_to_add:
                self.logger.warning('    %s' % seq_id)

        if failed_to_remove:
            self.logger.warning('  [Warning] Failed to remove the following sequence(s):')
            for seq_id in failed_to_remove:
                self.logger.warning('    %s' % seq_id)

        self.logger.info('')
        self.logger.info('  Modified genome written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - call_genes] Identifying genes within genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        prodigal = Prodigal(options.cpus, options.genes, options.output_dir)
        prodigal.run(genome_files)

        self.logger.info('')
        self.logger.info('  Identified genes written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - unique] Ensuring sequences are assigned to a single genome.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)

        duplicates = genome_tk.unique(genome_files)

        if len(duplicates) == 0:
            self.logger.info('  Pass: All sequences were identified exactly once.')
        else:
            self.logger.info('  Fail: One or more sequences were observed multiple times.')

            genome_ids = sorted(duplicates.keys())
            for i in xrange(0, len(genome_ids)):
                genome_idA = genome_ids[i]

                for j in xrange(i, len(genome_ids)):
                    genome_idB = genome_ids[j]

                    dup_seq_ids = duplicates[genome_idA][genome_idB]
                    if len(dup_seq_ids) == 0:
                        continue

                    self.logger.info('')
                    if genome_idA == genome_idB:
                        self.logger.info('  There are %d sequences present more than once in %s:' % (len(dup_seq_ids), genome_idA))
                    else:
                        self.logger.info('  There are %d sequences shared between %s and %s:' % (len(dup_seq_ids), genome_idA, genome_idB))

                    for seq_id in dup_seq_ids:
                        self.logger.info('    %s' % seq_id)

        self.time_keeper.print_time_stamp()

    def bin_compare(self, options):
        """Bin compare command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - bin_compare] Comparing two sets of genomes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.genomes_dir1)
        check_dir_exists(options.genomes_dir2)

        genomes_files1 = self._genome_files(options.genomes_dir1, options.genome_ext1)
        genomes_files2 = self._genome_files(options.genomes_dir2, options.genome_ext2)

        bin_comparer = BinComparer()
        bin_comparer.run(genomes_files1, genomes_files2, options.scaffold_file, options.output_file)

        self.logger.info('')
        self.logger.info('  Detailed bin comparison written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()

    def reference(self, options):
        """Reference command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - reference] Identifying scaffolds similar to specific genome(s).')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        # call genes on scaffold file
        prodigal = Prodigal(options.cpus)
        scaffold_genes_dir = os.path.join(options.output_dir, 'scaffold_genes')
        make_sure_path_exists(scaffold_genes_dir)
        prodigal.run([options.scaffold_file], False, 11, True, scaffold_genes_dir)
        scaffold_gene_file = os.path.join(scaffold_genes_dir, remove_extension(options.scaffold_file) + '.genes.faa')

        # call genes on reference genomes
        self.logger.info('')
        genome_files = self._genome_files(options.ref_genome_dir, options.genome_ext)
        ref_genome_genes_dir = os.path.join(options.output_dir, 'ref_genome_genes')
        make_sure_path_exists(ref_genome_genes_dir)
        prodigal.run(genome_files, False, None, False, ref_genome_genes_dir)

        ref_genome_gene_files = []
        for f in genome_files:
            ref_genome_gene_files.append(os.path.join(ref_genome_genes_dir, remove_extension(f) + '.genes.faa'))

        # get binned genomes
        binned_genomes = self._genome_files(options.ref_genome_dir, options.genome_ext)
        scaffold_id_bin_id = {}
        for g in binned_genomes:
            bin_id = remove_extension(g)
            for seq_id, _seq in seq_io.read_seq(g):
                #*** This is a bit questionable, but does make RefineM compatible with FinishM
                for contig_id in seq_id.split(':'):
                    scaffold_id_bin_id[contig_id] = bin_id
                scaffold_id_bin_id[seq_id] = bin_id

        reference = Reference(options.cpus, options.output_dir)
        reference_out = reference.run(options.scaffold_file,
                                        scaffold_gene_file,
                                        ref_genome_gene_files,
                                        options.db_file, options.coverage,
                                        options.evalue, options.per_identity, scaffold_id_bin_id)

        self.logger.info('')
        self.logger.info('  Results written to: ' + reference_out)

        self.time_keeper.print_time_stamp()

    def unbinned(self, options):
        """Unbinned Command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - unbinned] Identify unbinned scaffolds.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        check_dir_exists(options.genome_dir)

        genomes_files = self._genome_files(options.genome_dir, options.genome_ext)

        unbinned = Unbinned()
        unbinned_seqs = unbinned.run(genomes_files, options.scaffold_file, options.min_seq_len)

        seq_io.write_fasta(unbinned_seqs, options.output_file)

        self.logger.info('')
        self.logger.info('  Unbinned scaffolds written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()

    def coverage(self, options):
        """Coverage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - coverage] Calculating coverage of scaffolds.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        genome_files = self._genome_files(options.genome_dir, options.genome_ext)

        coverage = Coverage(options.cpus)
        coverage.run(genome_files, options.bam_files, options.output_file, options.all_reads, options.min_align, options.max_edit_dist)

        self.logger.info('')
        self.logger.info('  Coverage information written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()

    def tetra(self, options):
        """Tetra usage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - tetra] Calculating tetranucleotide signature of scaffolds.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        # determine sequences in each genome
        genome_files = self._genome_files(options.genome_dir, options.genome_ext)

        bin_id_scaffold_id = {}
        for genome_file in genome_files:
            bin_id = remove_extension(genome_file)
            for seq_id, _seq in seq_io.read_seq(genome_file):
                bin_id_scaffold_id[seq_id] = bin_id

        # calculate tetranucleotide signatures
        tetra = Tetranucleotide(options.cpus)
        seq_signatures = tetra.run(options.scaffold_file)

        # write results
        fout = open(options.output_file, 'w')

        fout.write('Scaffold id\tBin id')
        for kmer in tetra.canonical_order():
            fout.write('\t' + kmer)
        fout.write('\n')

        for seq_id, tetra_signature in seq_signatures.iteritems():
            fout.write(seq_id + '\t' + bin_id_scaffold_id.get(seq_id, 'unbinned'))
            fout.write('\t' + '\t'.join(map(str, tetra_signature)))
            fout.write('\n')

        fout.close()

        self.logger.info('')
        self.logger.info('  Tetranucleotide signatures written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""
        try:
            if options.bVerbose:
                logging.basicConfig(format='', level=logging.DEBUG)
            elif options.bQuiet:
                logging.basicConfig(format='', level=logging.ERROR)
            else:
                logging.basicConfig(format='', level=logging.INFO)
        except:
            logging.basicConfig(format='', level=logging.INFO)

        check_dependencies(('prodigal', 'diamond', 'ktImportText'))

        if(options.subparser_name == 'taxa_profile'):
            self.taxa_profile(options)
        elif(options.subparser_name == 'gene_profile'):
            self.gene_profile(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'bin_compare'):
            self.bin_compare(options)
        elif(options.subparser_name == 'reference'):
            self.reference(options)
        elif(options.subparser_name == 'call_genes'):
            self.call_genes(options)
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'tetra'):
            self.tetra(options)
        else:
            self.logger.error('  [Error] Unknown RefineM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
