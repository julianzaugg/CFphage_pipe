#!/usr/bin/env python
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
__author__ = "Julian Zaugg"
__copyright__ = "Copyright 2022"
__credits__ = ["Julian Zaugg"]
__license__ = "GPL3"
__maintainer__ = "Julian Zaugg"
__email__ = "j.zaugg near uq.edu.au"
__status__ = "Development"

###############################################################################
import argparse
import subprocess
import os
import sys
import signal
import pandas as pd
from yaml import dump as yaml_dump

from cfphage_pipe.__init__ import __version__


def phelp():
    print(
        """
                            ......:::::: CFPhage_pipe ::::::......
        Pipeline for the processing of (Cystic Fibrosis) isolates for phage discovery and evaluation
                assemble         - Assemble short, short + long, or just long reads for one or more isolates
                predict_virus    - Predict viruses in provided contigs
                annotate_isolate - Functionally annotate one or more isolates
                annotate_virus   
        """
    )


def str2bool(v):
    if isinstance(v, bool):
        return (v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return (True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return (False)
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# def get_snakefile_path(name):
#     thisdir = os.path.dirname(__file__)
#     snakefile = os.path.join(thisdir, "conf", name)
#     return snakefile


##########################################################################################
# Functions to source the conda environment variables
# Taken from https://github.com/rhysnewell/aviary/blob/93f7e342e2659562d8ae7a40816e46b47cea01e8/aviary/config/config.py
"""
Load the reference package. This will fail if the directory doesn't exist.
"""


def get_software_db_path(db_name='CONDA_ENV_PATH', software_flag='--conda-prefix'):
    try:
        source_conda_env()
        SW_PATH = os.environ[db_name]
        return SW_PATH
    except KeyError:
        try:
            source_conda_env()
            CONDA_PATH = os.environ[db_name]
            return CONDA_PATH
        except KeyError:
            try:
                source_bashrc()
                CONDA_PATH = os.environ[db_name]
                return CONDA_PATH
            except KeyError:
                print('\n' + '=' * 100)
                print(' ERROR '.center(100))
                print('_' * 100 + '\n')
                print(f"The '{db_name}' environment variable is not defined.".center(100) + '\n')
                print('Please set this variable to your default server/home directory conda environment path.'.center(
                    100))
                print(f'Alternatively, use {software_flag} flag.'.center(100))
                print('=' * 100)
                signal.alarm(120)
                os.environ[db_name] = input(f'Input {db_name} now:')
                try:
                    subprocess.Popen(
                        'mkdir -p %s/etc/conda/activate.d/; mkdir -p %s/etc/conda/deactivate.d/; echo "export %s=%s" >> %s/etc/conda/activate.d/cfphage_pipe.sh; echo "unset %s" >> %s/etc/conda/deactivate.d/cfphage_pipe.sh; ' %
                        (os.environ['CONDA_PREFIX'], os.environ['CONDA_PREFIX'], db_name, os.environ[db_name],
                         os.environ['CONDA_PREFIX'], db_name, os.environ['CONDA_PREFIX']), shell=True).wait()
                except KeyError:
                    subprocess.Popen(
                        'echo "export %s=%s" >> ~/.bashrc' %
                        (db_name, os.environ[db_name]), shell=True).wait()
                signal.alarm(0)
                print('=' * 100)
                print(
                    'Reactivate your cfphage_pipe conda environment or source ~/.bashrc to suppress this message.'.center(
                        100))
                print('=' * 100)

                return os.environ[db_name]


def source_conda_env():
    try:
        with open(format('%s/etc/conda/activate.d/cfpage_pipe.sh' % os.environ['CONDA_PREFIX'])) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                # if 'export' not in line:
                #     continue
                # Remove leading `export `, if you have those
                # then, split name / value pair
                # key, value = line.replace('export ', '', 1).strip().split('=', 1)
                try:
                    key, value = line.strip().split('=', 1)
                    os.environ[key] = value  # Load to local environ
                except ValueError:
                    continue
    except FileNotFoundError:
        # File not found so going to have to create it
        pass


def source_bashrc():
    try:
        with open('%s/.bashrc' % os.environ['HOME']) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                # if 'export' not in line:
                #     continue
                # Remove leading `export `, if you have those
                # then, split name / value pair
                # key, value = line.replace('export ', '', 1).strip().split('=', 1)
                try:
                    key, value = line.strip().split('=', 1)
                    os.environ[key] = value  # Load to local environ
                except ValueError:
                    continue
    except FileNotFoundError:
        # File not found so going to have to create it
        pass


##########################################################################################


def create_config(configfile, args):
    # Check databases
    if not os.path.exists(args.gtdbtk_db):
        print(f"Error: path to GTDB-Tk database {args.gtdbtk_db} does not exits")
        sys.exit(1)

    if not os.path.exists(args.imgvr_diamond_db):
        print(f"Error: path to IMGVR protein DIAMOND database {args.imgvr_db} does not exits")
        sys.exit(1)

    if not os.path.exists(args.amrfinderplus_db):
        print(f"Error: path to AMRFinderPlus database {args.amrfinderplus_db} does not exits")
        sys.exit(1)

    if not os.path.exists(args.checkv_db):
        print(f"Error: path to CheckV database {args.checkv_db} does not exits")
        sys.exit(1)

    if not os.path.exists(args.vibrant_db):
        print(f"Error: path to VIBRANT database {args.vibrant_db} does not exits")
        sys.exit(1)


def main():
    ############################### Main parser ###############################

    main_parser = argparse.ArgumentParser(prog='cfphage_pipe',
                                          formatter_class=CustomHelpFormatter,
                                          add_help=False)
    main_parser.add_argument('--version',
                             action='version',
                             version=__version__,
                             help='Show version information.')

    # snakemake.snakemake(snakefile, configfiles=[args.configfile])
    # parser=argparse.ArgumentParser(description="", usage="")
    # parser.add_argument('-c', '--configfile', action='')
    # args = parser.parse_args()
    ############################## Command groups #############################
    base_group = argparse.ArgumentParser(formatter_class=CustomHelpFormatter, add_help=False)

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
        default=100,
    )

    base_group.add_argument(
        '-o', '--output',
        help='Output directory',
        dest='output',
        default='./',
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
        '--build',
        help='Build conda environments and then exits. Equivalent to \"--snakemake-cmds \'--conda-create-envs-only True \' \"',
        type=str2bool,
        nargs='?',
        const=True,
        dest='build',
    )

    base_group.add_argument(
        '--conda-prefix', '--conda_prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default=get_software_db_path('CONDA_ENV_PATH', '--conda-prefix'),
    )

    base_group.add_argument(
        '--snakemake-cmds',
        help='Additional commands to supplied to snakemake in the form of a single string'
             'e.g. "--print-compilation True". '
             'NOTE: Most commands in snakemake -h are valid but some commands may clash with commands '
             'cfphage directly supplies to snakemake. Please make'
             "sure your additional commands don't clash.",
        dest='cmds',
        default='',
    )

    ############################### Parsing input ###############################

    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    else:
        args = main_parser.parse_args()

    prefix = args.output
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    if args.build:
        try:
            args.cmds = args.cmds + '--conda-create-envs-only '
        except TypeError:
            args.cmds = '--conda-create-envs-only '


############################### Classes ###############################
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
