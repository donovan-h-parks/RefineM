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
import traceback
from collections import defaultdict

import pysam

from biolib.common import remove_extension

from refinem.errors import ParsingError


class ReadLoader:
    """Callback for counting aligned reads with pysam.fetch"""

    def __init__(self, all_reads, min_align_per, max_edit_dist_per):
        self.all_reads = all_reads
        self.min_align_per = min_align_per
        self.max_edit_dist_per = max_edit_dist_per

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
        elif read.alen < self.min_align_per * read.rlen:
            self.num_failed_align_len += 1
        elif read.opt('NM') > self.max_edit_dist_per * read.rlen:
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
        self.logger = logging.getLogger('timestamp')
        self.reporter = logging.getLogger('no_timestamp')

        self.cpus = cpus

    def run(self, bam_files, out_file, all_reads, min_align_per, max_edit_dist_per):
        """Calculate coverage of sequences for each BAM file."""

        # make sure all BAM files are indexed
        for bam_file in bam_files:
            if not os.path.exists(bam_file + '.bai'):
                self.logger.error('BAM index file is missing: ' + bam_file + '.bai\n')
                sys.exit()

        # calculate coverage of each BAM file
        coverage_info = {}
        for i, bam_file in enumerate(bam_files):
            self.logger.info('Calculating coverage profile for %s (%d of %d):' % (ntpath.basename(bam_file), i + 1, len(bam_files)))

            coverage_info[bam_file] = mp.Manager().dict()
            coverage_info[bam_file] = self._process_bam(bam_file, all_reads, min_align_per, max_edit_dist_per, coverage_info[bam_file])

        fout = open(out_file, 'w')
        header = 'Scaffold Id\tLength (bp)'
        for bam_file in bam_files:
            bam_id = remove_extension(bam_file)
            header += '\t' + bam_id
        fout.write(header + '\n')

        for seq_id in coverage_info[coverage_info.keys()[0]].keys():
            row_str = seq_id + '\t' + str(coverage_info[coverage_info.keys()[0]][seq_id].seq_len)
            for bam_file in bam_files:
                bam_id = remove_extension(bam_file)
                if seq_id in coverage_info[bam_file]:
                    row_str += '\t' + str(coverage_info[bam_file][seq_id].coverage)
                else:
                    row_str += '\t' + '0.0'
            fout.write(row_str + '\n')

        fout.close()

    def _process_bam(self, bam_file, all_reads, min_align_per, max_edit_dist_per, coverage_info):
        """Calculate coverage of scaffolds in BAM file."""

        # determine coverage for each reference scaffolds
        worker_queue = mp.Queue()
        writer_queue = mp.Queue()

        bamfile = pysam.Samfile(bam_file, 'rb')
        ref_seq_ids = bamfile.references
        ref_seq_lens = bamfile.lengths

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
            worker_proc = [mp.Process(target=self._worker, args=(bam_file, all_reads, min_align_per, max_edit_dist_per, worker_queue, writer_queue)) for _ in range(self.cpus)]
            write_proc = mp.Process(target=self._writer, args=(coverage_info, len(ref_seq_ids), writer_queue))

            write_proc.start()

            for p in worker_proc:
                p.start()

            for p in worker_proc:
                p.join()

            writer_queue.put((None, None, None, None, None, None, None, None, None, None, None))
            write_proc.join()
        except:
            print traceback.format_exc()
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
                num_reads = 0
                num_mapped_reads = 0
                num_duplicates = 0
                num_secondary = 0
                num_failed_qc = 0
                num_failed_align_len = 0
                num_failed_edit_dist = 0
                num_failed_proper_pair = 0
                coverage = 0

                for read in bamfile.fetch(seq_id, 0, seq_len):
                    num_reads += 1

                    if read.is_unmapped:
                        pass
                    elif read.is_duplicate:
                        num_duplicates += 1
                    elif read.is_secondary or read.is_supplementary:
                        num_secondary += 1
                    elif read.is_qcfail:
                        num_failed_qc += 1
                    elif read.query_alignment_length < min_align_per * read.query_length:
                        num_failed_align_len += 1
                    elif read.get_tag('NM') > max_edit_dist_per * read.query_length:
                        num_failed_edit_dist += 1
                    elif not all_reads and not read.is_proper_pair:
                        num_failed_proper_pair += 1
                    else:
                        num_mapped_reads += 1

                        # Note: the alignment length (query_alignment_length) is used instead of the
                        # read length (query_length) as this bring the calculated coverage
                        # in line with 'samtools depth' (at least when the min
                        # alignment length and edit distance thresholds are zero)
                        coverage += read.query_alignment_length

                coverage = float(coverage) / seq_len

                queue_out.put((seq_id, seq_len, coverage, num_reads,
                                num_duplicates, num_secondary, num_failed_qc,
                                num_failed_align_len, num_failed_edit_dist,
                                num_failed_proper_pair, num_mapped_reads))

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

            if not self.logger.is_silent:
                processed_ref_seqs += 1
                statusStr = '  Finished processing %d of %d (%.2f%%) reference sequences.' % (processed_ref_seqs, num_reference_seqs, float(processed_ref_seqs) * 100 / num_reference_seqs)
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

            coverage_info[seq_id] = CoverageStruct(seq_len=seq_len, mapped_reads=num_mapped_reads, coverage=coverage)

        if not self.logger.is_silent:
            sys.stderr.write('\n')

        self.reporter.info('')
        self.reporter.info('  # total reads: %d' % total_reads)
        self.reporter.info('    # properly mapped reads: %d (%.1f%%)' % (total_mapped_reads, float(total_mapped_reads) * 100 / total_reads))
        self.reporter.info('    # duplicate reads: %d (%.1f%%)' % (total_duplicates, float(total_duplicates) * 100 / total_reads))
        self.reporter.info('    # secondary reads: %d (%.1f%%)' % (total_secondary, float(total_secondary) * 100 / total_reads))
        self.reporter.info('    # reads failing QC: %d (%.1f%%)' % (total_failed_qc, float(total_failed_qc) * 100 / total_reads))
        self.reporter.info('    # reads failing alignment length: %d (%.1f%%)' % (total_failed_align_len, float(total_failed_align_len) * 100 / total_reads))
        self.reporter.info('    # reads failing edit distance: %d (%.1f%%)' % (total_failed_edit_dist, float(total_failed_edit_dist) * 100 / total_reads))
        self.reporter.info('    # reads not properly paired: %d (%.1f%%)' % (total_failed_proper_pair, float(total_failed_proper_pair) * 100 / total_reads))

    def read(self, coverage_file):
        """Read coverage information from file.

        Parameters
        ----------
        coverage_file : str
            File containing coverage profiles.

        Returns
        -------
        dict : d[scaffold_id][bam_id] -> coverage
            Coverage profile for each scaffold.
        dict : d[scaffold_id] -> length
            Length of each scaffold.
        """

        try:
            coverage = defaultdict(lambda: defaultdict(float))
            length = {}
            with open(coverage_file) as f:
                header = f.readline().split('\t')
                bam_ids = [x.strip() for x in header[2:]]

                for line in f:
                    line_split = line.split('\t')
                    scaffold_id = line_split[0]
                    scaffold_len = int(line_split[1])

                    length[scaffold_id] = scaffold_len

                    for i, cov in enumerate(line_split[2:]):
                        coverage[scaffold_id][bam_ids[i]] = float(cov)
        except IOError:
            self.logger.error('Failed to open signature file: %s' % coverage_file)
            sys.exit()
        except:
            print traceback.format_exc()
            print ''
            raise ParsingError("[Error] Failed to process coverage file: " + coverage_file)
            sys.exit()

        return coverage, length
