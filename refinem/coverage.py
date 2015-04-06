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

import sys
import os
import multiprocessing as mp
import logging
import ntpath
from collections import defaultdict

import pysam

import biolib.seq_io as seq_io
from biolib.common import remove_extension

from numpy import mean, std


class ReadLoader:
    """Callback for counting aligned reads with pysam.fetch"""

    def __init__(self, all_reads, min_align_per, max_edit_dist_per):
        self.all_reads = all_reads
        self.min_align_per = min_align_per
        self.max_edit_dist_per = max_edit_dist_per

        self.num_reads = 0
        self.num_mapped_reads = 0
        self.num_duplicates = 0
        self.num_secondary = 0
        self.num_failed_qc = 0
        self.num_failed_align_len = 0
        self.num_failed_edit_dist = 0
        self.num_failed_proper_pair = 0

        self.coverage = 0

    def __call__(self, read):
        self.num_reads += 1

        if read.is_unmapped:
            pass
        elif read.is_duplicate:
            self.num_duplicates += 1
        elif read.is_secondary:
            self.num_secondary += 1
        elif read.is_qcfail:
            self.num_failed_qc += 1
        elif read.alen < self.min_align_per*read.rlen:
            self.num_failed_align_len += 1
        elif read.opt('NM') > self.max_edit_dist_per*read.rlen:
            self.num_failed_edit_dist += 1
        elif not self.all_reads and not read.is_proper_pair:
            self.num_failed_proper_pair += 1
        else:
            self.num_mapped_reads += 1

            # Note: the alignment length (alen) is used instead of the
            # read length (rlen) as this bring the calculated coverage
            # in line with 'samtools depth' (at least when the min
            # alignment length and edit distance thresholds are zero).
            self.coverage += read.alen

class CoverageStruct():
    """Coverage information for scaffolds."""
    
    def __init__(self, seq_len, mapped_reads, coverage):
        self.seq_len = seq_len
        self.mapped_reads = mapped_reads
        self.coverage = coverage

class Coverage():
    """Calculate coverage of all sequences."""
    
    def __init__(self, cpus):
        self.logger = logging.getLogger()
        
        self.cpus = cpus 

    def run(self, genome_files, bam_files, out_file, all_reads, min_align_per, max_edit_dist_per):
        """Calculate coverage of sequences for each BAM file."""
        
        # determine genome assignment of each scaffold
        self.logger.info('  Determining genome assignment of each scaffold.')
            
        seq_id_genome_id = {}
        for genome_file in genome_files:
            genome_id = remove_extension(genome_file)
            for seq_id, _seq in read_seq(genome_file):
                seq_id_genome_id[seq_id] = genome_id
        
        # process each fasta file
        self.logger.info("  Processing %d file(s) with %d cpus.\n" % (len(bam_files), self.cpus))
            
        # make sure all BAM files are sorted
        for bam_file in bam_files: 
            if not os.path.exists(bam_file + '.bai'):
                self.logger.error('  [Error] BAM file is not sorted: ' + bam_file + '\n')
                sys.exit()
 
        # calculate coverage of each BAM file
        coverage_info = {}
        for i, bam_file in enumerate(bam_files): 
            self.logger.info('  Processing %s (%d of %d):' % (ntpath.basename(bam_file), i+1, len(bam_files)))
            
            coverage_info[bam_file] = mp.Manager().dict()
            coverage_info[bam_file] = self._process_bam(bam_file, all_reads, min_align_per, max_edit_dist_per, coverage_info[bam_file]) 
   
        fout = open(out_file, 'w')
        header = 'Scaffold Id\tGenome Id\tScaffold length (bp)'
        for _ in bam_files:
            header += '\tBam Id\tCoverage\tMapped reads'
        
        fout.write(header + '\n')

        for seq_id in coverage_info[coverage_info.keys()[0]].keys():
            row_str = seq_id + '\t' + seq_id_genome_id.get(seq_id, 'unbinned') + '\t' + str(coverage_info[coverage_info.keys()[0]][seq_id].seq_len)
            for bam_file in bam_files:
                bam_id = remove_extension(bam_file)
                row_str += '\t' + bam_id + '\t' + str(coverage_info[bam_file][seq_id].coverage) + '\t' + str(coverage_info[bam_file][seq_id].mapped_reads)
            fout.write(row_str + '\n')
            
        fout.close()

    def _process_bam(self, bam_file, all_reads, min_align_per, max_edit_dist_per, coverage_info):
        """Calculate coverage of scaffolds in BAM file."""

        # determine coverage for each reference scaffolds
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        bam_file = pysam.Samfile(bam_file, 'rb')
        ref_seq_ids = bam_file.references
        ref_seq_lens = bam_file.lengths

        # populate each thread with reference scaffolds to process
        # Note: reference scaffolds are sorted by number of mapped reads
        # so it is important to distribute reads in a sensible way to each
        # of the threads
        ref_seq_lists = [[] for _ in range(self.cpus)]
        ref_len_lists = [[] for _ in range(self.cpus)]

        cpu_index = 0
        incDir = 1
        for ref_seq_id, ref_len in zip(ref_seq_ids, ref_seq_lens):
            ref_seq_lists[cpu_index].append(ref_seq_id)
            ref_len_lists[cpu_index].append(ref_len)

            cpu_index += incDir
            if cpu_index == self.cpus:
                cpu_index = self.cpus - 1
                incDir = -1
            elif cpu_index == -1:
                cpu_index = 0
                incDir = 1

        for i in range(self.cpus):
            worker_queue.put((ref_seq_lists[i], ref_len_lists[i]))

        for _ in range(self.cpus):
            worker_queue.put((None, None))

        try:
            worker_proc = [mp.Process(target = self._worker, args = (bam_file, all_reads, min_align_per, max_edit_dist_per, worker_queue, writer_queue)) for _ in range(self.cpus)]
            write_proc = mp.Process(target = self._writer, args = (coverage_info, len(ref_seq_ids), writer_queue))

            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put((None, None, None, None, None, None, None, None, None, None, None))
            write_proc.join()
        except:
            for p in worker_proc:
                p.terminate()

            write_proc.terminate()
               
        return coverage_info
                   
    def _worker(self, bam_file, all_reads, min_align_per, max_edit_dist_per, queue_in, queue_out):
        """Process scaffold in parallel.
        
        Parameters
        ----------
        bam_file : str
            BAM file to process.
        all_reads : boolean
            Flag indicating if all reads or just paired reads should be processed.
        min_align_per : float
            Alignment percentage threshold for accepting mapped reads.
        max_edit_dist_per : float
            Edit distance threshold for accepting mapped reads.
        queue_in : queue
            Queue containing reference sequences and their lengths.
        queue_out : queue
            Queue to hold coverage results.
        """
        while True:
            seq_ids, seq_lens = queue_in.get(block=True, timeout=None)
            if seq_ids == None:
                break

            bamfile = pysam.Samfile(bam_file, 'rb')

            for seq_id, seq_len in zip(seq_ids, seq_lens):
                readLoader = ReadLoader(seq_len, all_reads, min_align_per, max_edit_dist_per)
                bamfile.fetch(seq_id, 0, seq_len, callback = readLoader)

                coverage = float(readLoader.coverage) / seq_len

                queue_out.put((seq_id, seq_len, coverage, readLoader.num_reads, 
                                readLoader.num_duplicates, readLoader.num_secondary, readLoader.num_failed_qc, 
                                readLoader.num_failed_align_len, readLoader.num_failed_edit_dist, 
                                readLoader.num_failed_proper_pair, readLoader.num_mapped_reads))

            bamfile.close()

    def _writer(self, coverage_info, num_reference_seqs, writer_queue):
        """Record coverage information for each scaffold.
        
        Parameters
        ----------
        coverage_info : managed dictionary
            Dictionary for recording CoverageStruct statistics for each scaffold.
        num_reference_seqs : int
            Number of reference scaffolds to process.
        writer_queue : queue
            Queue contain results of worker threads.
        """
        total_reads = 0
        total_duplicates = 0
        total_secondary = 0
        total_failed_qc = 0
        total_failed_align_len = 0
        total_failed_edit_dist = 0
        total_failed_proper_pair = 0
        total_mapped_reads = 0

        processed_ref_seqs = 0
        while True:
            seq_id, seq_len, coverage, num_reads, num_duplicates, num_secondary, num_failed_qc, num_failed_align_len, num_failed_edit_dist, num_failed_proper_pair, num_mapped_reads = writer_queue.get(block=True, timeout=None)
            if seq_id == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processed_ref_seqs += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) reference sequences.' % (processed_ref_seqs, num_reference_seqs, float(processed_ref_seqs)*100/num_reference_seqs)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

                total_reads += num_reads
                total_duplicates += num_duplicates
                total_secondary += num_secondary
                total_failed_qc += num_failed_qc
                total_failed_align_len += num_failed_align_len
                total_failed_edit_dist += num_failed_edit_dist
                total_failed_proper_pair += num_failed_proper_pair
                total_mapped_reads += num_mapped_reads
                
            coverage_info[seq_id] = CoverageStruct(seq_len = seq_len, mapped_reads = num_mapped_reads, coverage = coverage)

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')
            
            print ''
            print '    # total reads: %d' % total_reads
            print '      # properly mapped reads: %d (%.1f%%)' % (total_mapped_reads, float(total_mapped_reads)*100/total_reads)
            print '      # duplicate reads: %d (%.1f%%)' % (total_duplicates, float(total_duplicates)*100/total_reads)
            print '      # secondary reads: %d (%.1f%%)' % (total_secondary, float(total_secondary)*100/total_reads)
            print '      # reads failing QC: %d (%.1f%%)' % (total_failed_qc, float(total_failed_qc)*100/total_reads)
            print '      # reads failing alignment length: %d (%.1f%%)' % (total_failed_align_len, float(total_failed_align_len)*100/total_reads)
            print '      # reads failing edit distance: %d (%.1f%%)' % (total_failed_edit_dist, float(total_failed_edit_dist)*100/total_reads)
            print '      # reads not properly paired: %d (%.1f%%)' % (total_failed_proper_pair, float(total_failed_proper_pair)*100/total_reads)
            print ''
            
    def parse_coverage(self, coverage_file):
        """Read coverage information from file.
        
        Parameters
        ----------
        coverage_file : str
            File containing coverage profiles.
        
        Returns
        -------
        dict : d[genome_id][scaffold_id][bam_id] -> coverage
            Coverage profile for each genome.
        """
        
        coverage_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        with open(coverage_file) as f:
            f.readline()
            
            for line in f:
                line_split = line.split('\t')
                scaffold_id = line_split[0]
                genome_id = line_split[1]
  
                for i in xrange(3, len(line_split), 3):
                    bam_id = line_split[i]
                    coverage = float(line_split[i+1])
                    coverage_stats[genome_id][scaffold_id][bam_id] = coverage
                
        return coverage_stats
                    