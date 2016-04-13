#!/usr/bin/env python

import pysam, sys, os, glob
from jobTree.src.bioio import reverseComplement, fastaRead, system, fastaWrite, \
cigarRead, logger, nameValue
import numpy as np

# file 1 is meta seq file
# make scratch dir
system("cPecanMultipleAlign %s" %(sys.argv[1])) 



