"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import time
import argparse
import os
import pickle

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate the global trajectory object file from the cdf file.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        required=True,
                        help='[Required] The 3-letters code of the ligand',
                        nargs=1)
    parser.add_argument('-n',
                        dest='gt_name',
                        help='[Optional] The name of the gt object file (Default: name of the cdf file)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the gt object file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-alh',
                        dest='alh',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Generate atom-level hydrophobic features')
    
    return parser.parse_args()

if __name__ == '__main__':
    args = parseArguments()

    initial_time = time.time()
    cdf = args.cdf[0]
    ligand_code = args.ligand_code[0]

    if args.gt_name is None:
        gt_name = os.path.basename(cdf).split('_chunk_')[0][:-4]
    else:
        gt_name = args.gt_name[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    ph4_interaction_dictionary = getPh4InteractionDictionary(cdf, ligand_code, args.alh)

    with open(output + gt_name + '.gt', 'wb') as handle:
        pickle.dump(ph4_interaction_dictionary, handle)

    calc_time = time.time() - initial_time
    print('> Global trajectory object file generated in {}s'.format(int(calc_time)))
