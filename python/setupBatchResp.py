#!/usr/bin/python

#=============================================================================================
# Program to test the Resp.py wrapper for the RESP charge derivation procedure.
#=============================================================================================
# TODO:
# - add command-line arguments
#=============================================================================================
# Version control information
__version__ = "$Revision: 1.1 $"
# $Source: /share_kubo/Storage/cvsroot/pymol-tools/setupBatchResp.py,v $
#=============================================================================================
# Imports
import sys
import os.path
from optparse import OptionParser # For parsing of command line arguments
from openeye.oechem import *       # Use OpenEye OEChem toolkit
from openeye.oeomega import *      # Use OpenEye Omega
import BatchResp # for RESP charge-fitting
import exceptions
#=============================================================================================

# Create command-line argument options.
usage_string = """\
usage: %prog --input INFILE.mol2 --outdir OUTPUT_DIRECTORY

example: %prog --input guthrie.mol2 --outdir molecules/
"""

version_string = "%prog %__version__"

parser = OptionParser(usage=usage_string, version=version_string)

parser.add_option("-i", "--input", metavar='INFILE',
                  action="store", type="string", dest='infile', default='',
                  help="Input mol2 file containing one or more molecules.")
parser.add_option("-o", "--outdir", metavar='OUTPUT_DIRECTORY',
                  action="store", type="string", dest='output_directory', default='.',
                  help="Directory to contain output mol2 files, whose names are extracted from input mol2 file (default '.').")

# Parse command-line arguments.
(options,args) = parser.parse_args()

# Perform minimal error checking.
if not options.infile:
  parser.print_help()
  parser.error("An input file must be specified.")

# Open input mol2 file.
input_molecule_stream = oemolistream()
input_molecule_stream.open(options.infile)

# Create a molecule.
molecule = OEMol()

# Create an instance of Resp, telling it where to find GAMESS and Antechamber.
#resp = Resp(gamess_path = '/Users/jchodera/local/src/gamess', antechamber_path = '/Users/jchodera/local/src/antechamber-1.27')
batch_resp = BatchResp.BatchResp()

# Keep track of number of conformers for each.
molecules = []

# Process all molecules in the input stream
while OEReadMolecule(input_molecule_stream, molecule):
    # Extract name of molecule.
    molecule_name = molecule.GetTitle().split()[0].replace('.xyz','')

    print "\nProcessing '%s'..." % molecule_name

    # Make a directory for this molecule.
    molecule_directory = os.path.join(options.output_directory, molecule_name)
    print "Creating directory '%s'..." % molecule_directory
    os.mkdir(molecule_directory)
    
    # Set up GAMESS directories.
    nconformers = batch_resp.setupJobs(molecule, molecule_name, molecule_directory, generate_conformers = True)

    # Keep track of number of conformers.
    molecules.append( (molecule_name, nconformers) )

# Generate master submission script in order of fewest conformers to most conformers.
outfile = open(os.path.join(options.output_directory,'submit_all.sh'),'w')
import operator
molecules.sort(key = operator.itemgetter(1))
for molecule in molecules:
    print molecule
    outfile.write("cd %s ; qsub run_batch.cmd ; cd .. \n" % molecule[0])
    print "%5d %s" % (molecule[1], molecule[0])
outfile.close()


    
    
    
