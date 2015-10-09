#!/usr/bin/python

#=============================================================================================
# Program to analyze a completed batch of quantum chemical calculations set up by setupBatchResp.py.
#=============================================================================================
# TODO:
#=============================================================================================
# Version control information
__version__ = "$Revision: 1.2 $"
# $Source: /share_kubo/Storage/cvsroot/pymol-tools/analyzeBatchResp.py,v $
#=============================================================================================
# Imports
import sys
import os
from optparse import OptionParser # For parsing of command line arguments
from openeye.oechem import *       # Use OpenEye OEChem toolkit
from openeye.oeomega import *      # Use OpenEye Omega
import BatchResp # for RESP charge-fitting
import exceptions
#=============================================================================================

# Create command-line argument options.
usage_string = """\
usage: %prog --dir DIRECTORY

example: %prog --dir molecules/ch3oph
"""

version_string = "%prog %__version__"

parser = OptionParser(usage=usage_string, version=version_string)

parser.add_option("-d", "--dir", metavar='INFILE',
                  action="store", type="string", dest='dir', default='',
                  help="Directory in which QM calcs were performed.")

# Parse command-line arguments.
(options,args) = parser.parse_args()

# Perform minimal error checking.
if not options.dir:
  parser.print_help()
  parser.error("A directory must be specified")

# Create an instance of Resp, telling it where to find GAMESS and Antechamber.
batch_resp = BatchResp.BatchResp()

# Perform RESP fit.
batch_resp.respFit(options.dir, net_charge = 0)
    
