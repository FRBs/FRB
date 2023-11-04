""" Generate a set of fake FRBs for testing """

import numpy as np

from zdm.chime import grids

from IPython import embed

def for_chime():


    all_rates, all_singles, all_reps = grids.load()
    embed(header='14 of mock.py')



if __name__ == '__main__':
    for_chime()