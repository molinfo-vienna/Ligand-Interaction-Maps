"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import math
import sys
import time
import os
import copy

import MDAnalysis
import argparse
import numpy

import scipy.stats
import warnings
from collections import defaultdict
from mpl_toolkits.axes_grid1 import make_axes_locatable
from multiprocessing.dummy import Pool
import matplotlib
import matplotlib.pylab

import CDPL

if __name__ == '__main__':
    print('> All dependencies are satisfied.')
