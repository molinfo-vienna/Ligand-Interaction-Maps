"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import argparse
import os

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a cdf file (CDPKIT format) from both a topology and a trajectory file.')
    parser.add_argument('-trj',
                        dest='trajectory',
                        required=True,
                        help='[Required] The path of the trajectory file, can also handle multiple trajectories of the same system',
                        nargs="*")
    parser.add_argument('-top',
                        dest='topology',
                        required=True,
                        help='[Required] The path of the topology file',
                        nargs=1)
    parser.add_argument('-n',
                        dest='cdf_name',
                        help='[Optional] The name of the cdf file (Default: name of the dcd file)',
                        nargs=1,
                        default=None)
    parser.add_argument('-cs',
                        dest='chunk_size',
                        help='[Optional] The number of frames to consider per chunk (default: no chunks)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    return parser.parse_args()

if __name__ == '__main__':
    args = parseArguments()

    trajectory = args.trajectory
    topology = args.topology[0]

    if args.cdf_name is None:
        cdf_name = os.path.basename(trajectory[0])[:-4]
    else:
        cdf_name = args.cdf_name[0]

    if args.chunk_size is None:
        chunk_size = 500
    else:
        chunk_size = int(args.chunk_size[0])

    if args.output is None:
        output = './'
    else:
        output = args.output[0] + '/'

    if output[-1] != '/':
        output += '/'

    chunk_number = cdfMol(topology, trajectory, output, cdf_name, chunk_size=chunk_size)
    mergeCDFMolecule(output + cdf_name + '_chunk_0.cdf', chunk_number)
