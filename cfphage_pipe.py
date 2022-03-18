# cfphage_pipe.py - Info about cfphage_pipe.py
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
from aviary.__init__ import __version__
import aviary.config.config as Config
from aviary.modules.processor import Processor
__author__ = "Julian Zaugg"
__copyright__ = "Copyright 2022"
__credits__ = ["Julian Zaugg"]
__license__ = "GPL3"
__maintainer__ = "Julian Zaugg"
__email__ = "j.zaugg near uq.edu.au"
__status__ = "Development"

###############################################################################
# System imports
import sys
import argparse
import logging
import os
import re
import tempfile
import glob
from datetime import datetime


def str2bool(v):
    if isinstance(v, bool):
        return(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return(True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return(False)
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    ############################ ~ Main Parser ~ ##############################
    main_parser = argparse.ArgumentParser(prog='cfphage_pipe',
                                          formatter_class=CustomHelpFormatter,
                                          add_help=False)
    main_parser.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')
    subparsers = main_parser.add_subparsers(help="--", dest='subparser_name')


    ####################################################################

    base_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                         add_help=False)

    base_group.add_argument(
        '-t', '--max-threads', '--max_threads',
        help='Maximum number of threads given to any particular process',
        dest='max_threads',
        default=8,
    )

    base_group.add_argument(
        '-m', '--max-memory', '--max_memory',
        help='Maximum memory for available usage in Gigabytes',
        dest='max_memory',
        default=250,
    )

    base_group.add_argument(
        '-o', '--output',
        help='Output directory',
        dest='output',
        default='./',
    )

    base_group.add_argument(
        '--conda-prefix', '--conda_prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default=Config.get_software_db_path('CONDA_ENV_PATH', '--conda-prefix'),
    )

    base_group.add_argument(
        '--dry-run', '--dry_run', '--dryrun',
        help='Perform snakemake dry run, tests workflow order and conda environments',
        type=str2bool,
        nargs='?',
        const=True,
        dest='dryrun',
        default=False,
    )

    base_group.add_argument(
        '--conda-frontend', '--conda_frontend',
        help='Which conda frontend to use, mamba is faster but harder to debug. Switch this to conda '
             'If experiencing problems installing environments',
        dest='conda_frontend',
        default="mamba",
        choices=["conda", "mamba"],
    )

    base_group.add_argument(
        '--build',
        help='Build conda environments and then exits. Equivalent to \"--snakemake-cmds \'--conda-create-envs-only True \' \"',
        type=str2bool,
        nargs='?',
        const=True,
        dest='build',
    )

    base_group.add_argument(
        '--snakemake-cmds',
        help='Additional commands to supplied to snakemake in the form of a single string'
             'e.g. "--print-compilation True". '
             'NOTE: Most commands in snakemake -h are valid but some commands may clash with commands '
             'aviary directly supplies to snakemake. Please make'
             "sure your additional commands don't clash.",
        dest='cmds',
        default='',
    )

    ####################################################################
    long_read_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter,
                                              add_help=False)
    long_read_group.add_argument(
        '-l', '--longreads', '--long-reads', '--long_reads',
        help='A space separated list of interleaved read files for the binning process. NOTE: If performing assembly and '
             'multiple long read files are provided, then only the first file is used for assembly. ',
        dest='longreads',
        nargs='*',
        default="none"
    )

###############################################################################
################################ - Classes - ##################################

class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])


if __name__ == '__main__':
    sys.exit(main())
