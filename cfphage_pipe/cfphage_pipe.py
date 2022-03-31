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

__author__ = "Julian Zaugg"
__copyright__ = "Copyright 2022"
__credits__ = ["Julian Zaugg"]
__license__ = "GPL3"
__maintainer__ = "Julian Zaugg"
__email__ = "j.zaugg near uq.edu.au"
__status__ = "Development"

###############################################################################
import argparse
import os

import snakemake

def str2bool(v):
    if isinstance(v, bool):
        return(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return(True)
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return(False)
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_snakefile_path(name):
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, "conf", name)
    return snakefile

# def main(args):
    # snakefile=snakefile_path
    # snakemake.snakemake(snakefile, configfiles=[args.configfile])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description="", usage="")
    parser.add_argument('-c', '--configfile', action='')
    args = parser.parse_argse()
