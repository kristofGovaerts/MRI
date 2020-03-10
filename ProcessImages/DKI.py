#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Kristof
#
# Created:     22/09/2013
# Copyright:   (c) Kristof 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
from scipy.stats import kurtosis

def main():
    pass

def kurt(obs):
    num = np.sum((obs - np.mean(obs)) ** 4)/ len(obs)
    denom=np.var(obs)**2
    return num/denom

kurtosis([1,2,3,4,5])

dist=[random.gauss(10,2) for i in range(10000)]




if __name__ == '__main__':
    main()

