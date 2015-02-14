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

from cleanm.taxonomic_profile import TaxonomicProfile
from biolib.misc.time_keeper import TimeKeeper
from biolib.common import make_sure_path_exists


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()

    def taxa_profile(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CleanM - taxa_profile] Generating taxonomic profiles.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        make_sure_path_exists(options.output_dir)

        genome_files = []
        for f in os.listdir(options.genome_dir):
            if f.endswith(options.genome_ext):
                genome_files.append(os.path.join(options.genome_dir, f))

        if not genome_files:
            self.logger.warning('  [Warning] No genomes found. Check the --genome_ext flag used to identify genomes.')
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
        self.logger.info('  Taxonomic profiles written to: %s' % options.output_dir)

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

        if(options.subparser_name == 'taxa_profile'):
            self.taxa_profile(options)
        else:
            self.logger.error('  [Error] Unknown AAI command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
