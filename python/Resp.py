#!/usr/bin/python

#=============================================================================================
# Module for assigning RESP charges using GAMESS.
# Written by John D. Chodera <jchodera@gmail.com> 2007-01-31
#=============================================================================================
# TODO:
#=============================================================================================
# Version control information
__version__ = "$Revision: 1.1 $"
# $Source: /share_kubo/Storage/cvsroot/pymol-tools/Resp.py,v $
#=============================================================================================
# Imports
import sys
import os
import os.path
import tempfile                    # For managing temporary directories
from openeye.oechem import *       # Use OpenEye OEChem toolkit
from openeye.oeomega import *      # Use OpenEye Omega
import exceptions
#=============================================================================================

class Resp:
    """
    A class for constructing RESP-derived HF 6-31G* charges using GAMESS and Antechamber.
    
    """

    def __init__(self, gamess_path = None, antechamber_path = None):
        """
        Initialize the RESP module, ensuring that GAMESS and Antechamber can be found.
        """

        # Guess GAMESS path if not specified.
        if (gamess_path == None):
            gamess_path = os.environ['GMSPATH']
        # Check for the existence of a working GAMESS installation.
        if not os.path.exists(os.path.join(gamess_path,'rungms')):
            raise ParameterException("GAMESS not found on provided path (" + gamess_path + "). Please provide a valid path via 'gamess_path' optional argument or set GMSPATH environment variable.")
        # Store path.
        self.gamess_path = gamess_path
        
        # Guess Antechamber path if not specified.
        if (antechamber_path == None):
            antechamber_path = os.environ['ACHOME']
        # Check for the existence of a working Antechamber installation.
        if not os.path.exists(os.path.join(antechamber_path,'exe','antechamber')):
            raise ParameterException("Antechamber not found on provided path (" + antechamber_path + "). Please provide a valid path via 'antechamber_path' optional argument or set ACHOME environment variable.")
        # Store path.
        self.antechamber_path = antechamber_path

    def writeFile(self, filename, contents):
        """
        Write the specified text contents to a file.

        required arguments:
          filename - the filename to write to, created or overwritten
          contents - the string containing the contents of the file to write
        """

        # Open file for writing in text mode.
        outfile = open(filename, 'wt')

        # Write contents of file.
        outfile.write(contents)

        # Close the file.
        outfile.close()

    def readFile(self, filename):
        """
        Read the text contents of the specified file.

        required arguments:
          filename - the filename of the file whose contents are to be read

        return values:
          lines - the contents of the file, as a list of lines.
        """
        
        # Open the file for reading in text mode.
        infile = open(filename, 'rt')

        # Read contents of file to a list.
        lines = infile.readlines()
          
        # Return lines.
        return lines

    def extractLines(self, lines, first_line, last_line, retain_first_line = False, retain_last_line = False):
        """
        Extract the set of lines between a given set of lines.

        required arguments:
          lines - the set of lines to search through
          first_line, last_line - the pair of lines bracketing the section of interest

        optional arguments:
          

        return values:
          selected_lines - the block of lines of interest

        """

        # TODO

        return

    # TODO: Change optimize_geometry to True
    def computePotential(self, molecule, calculation = 'optimize', workdir = None):
        """
        Compute the electrostatic potential on concentric Connolly surfaces around the molecule using GAMESS.

        required arguments:
          molecule - the conformer to compute the electrostatic potential for

        optional arguments:
          calculation - one of 'optimize' (hours), 'singlepoint' (minutes), or 'test' (seconds) (defaults to 'optimize')
          workdir - working directory to set up GAMESS calculations in (default is to create a temporary directory and delete it after calculation)

        return values:
          potential - the (optimized) coordinates and potential 

        """

        # GAMESS header blocks.
        # GAMESS is instructed to carry out an RHF claculation with a 6-31G* basis set, after which
        # the electrostatic potential on concentric Connolly surfaces is evaluated.
        
        # Header for energy minimization.
        gamess_header['optimize'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=OPTIMIZE COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS $END
 $SCF    DIRSCF=.TRUE. $END
 $STATPT NSTEP=100 NPRT=-2 NPUN=-2 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
 """

        # Header for no optimization case.
        gamess_header['singlepoint'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=ENERGY COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS $END
 $SCF    DIRSCF=.TRUE. $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
 """

        # Header for test case, where initial HUCKEL guess is used instead for speed.
        # This should not be used in production.
        gamess_header['test'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=PROP COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS $END
 $SCF    DIRSCF=.TRUE. $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
 """

        # Check calculation option.
        if not gamess_header.has_key(calculation):
            raise ParameterError("Optional argument 'calculation' must be one of: " + gamess_header.keys())
        
        # Construct the GAMESS input file.
        # Select header for appropriate calculation.
        gamess_input = gamess_header[calculation]
        # Append data block containing Cartesian atomic coordinates.
        gamess_input += "$DATA\n"  # block header
        gamess_input += molecule.GetTitle() + "\n"  # molecule name        
        gamess_input += "C1\n" # symmetry point group        
        for atom in molecule.GetAtoms():  # atomic coordinates
            xyz = molecule.GetCoords(atom)
            gamess_input += "%10s %4.0f %13f %13f %13f\n" % (atom.GetName(), atom.GetAtomicNum(), xyz[0], xyz[1], xyz[2])
        gamess_input += "$END\n"  # terminate block

        # Create a temporary directory to contain GAMESS input and output if not specified by user.
        delete_workdir = False
        if (workdir == None):
            workdir = tempfile.mkdtemp()
            delete_workdir = True
        print "Setting up GAMESS calculation in '%s'..." % workdir

        # DEBUG
        # Write mol2 file of molecule in temporary directory.
        output_molecule_stream = oemolostream()
        output_molecule_stream.open(os.path.join(workdir, 'molecule.mol2'))
        OEWriteMolecule(output_molecule_stream, molecule)
        output_molecule_stream.close()

        # Write input file to temporary directory.
        self.writeFile(os.path.join(workdir, 'gamess.inp'), gamess_input)

        # Run GAMESS from the temporary directory.
        # TODO: Trap system command so that failure or abort will not cause us to end up in the wrong directory.
        olddir = os.getcwd()
        os.chdir(workdir)
        command = 'gms gamess.inp -l gamess.log'
        print command # DEBUG
        os.system(command)
        os.chdir(olddir)

        # Read .log and .dat file output.
        gamess_log_lines = self.readFile(os.path.join(workdir,'gamess.log'))
        gamess_dat_lines = self.readFile(os.path.join(workdir,'gamess.dat'))

        # Parse the .dat file to determine the final atomic coordinates and electrostatic potential.
        import re
        natoms_pattern = re.compile('TOTAL NUMBER OF ATOMS')
        npoints_pattern = re.compile('NUMBER OF POINTS SELECTED FOR FITTING')
        coords_pattern = re.compile('^atom')
                     
        mode = 'none'
        for line in gamess_dat_lines:
            if line.find('TOTAL NUMBER OF ATOMS')
            if mode == 'store_coords':
                # Store the atomic coordinates.
                
            elif mode == 'store_esp':
#                # Store the electrostatic potential
#                # Format is: IPT,X,Y,Z,ELPOTT
#                (index, x, y, z, elpott) = split(line)
#            if store_coords
#            
#            if line[0:19] == :
#                store_lines = True
#            if line[0:8] ==             
#
#        # Extract section containing electrostatic potential from .dat file.
#        esp_lines = extractLines(' ELECTRIC POTENTIAL', '--------')
#        
#        # Extract section containing optimized coordinates from .dat file.
#        coord_lines = extractLines('-------- START OF -MOLPLT- INPUT FILE ----------', '-------- END OF -MOLPLT- INPUT FILE ----------')
#        # Only keep the lines that will contain the atom locations
#        coord_lines = coord_lines[5:(5+natoms-1)]

#        if delete_workdir:
        # Clean up temporary directory.
#        for filename in os.listdir(workdir):
#            os.remove(os.path.join(workdir,filename))
#        os.rmdir(workdir)

        # return potential
        return gamess_dat_lines

    def computeCharges(self, molecule, generate_conformers = False):
        """
        Compute RESP-derived charges for the given OEchem molecule, which may contain multiple conformers.

        required arguments:
          molecule - an OEChem molecule

        optional arguments:
          generate_conformers - if True, multiple conformers will be generated by OpenEye's Omega (default False)

        returns:
          charged_molecule - the molecule passed as input with RESP charges assigned
        """

        print 'computeCharges: "%s"' % molecule.GetTitle()
#        if molecule.GetTitle()[0:9] != 'bencl.xyz':
#            return

        # Generate multiple conformers if requested.
        if (generate_conformers):
            # Initialize an instance of Omega.
            omega = OEOmega()
            # Generate multiple conformers from the given molecule.
            omega(molecule)            

        # Perform electronic structure calculations on each conformer.
        for conformer in molecule.GetConfs():
            # Compute optimized geometry and electrostatic potential.
            potential = self.computePotential(conformer, optimize_geometry = False)

            # TODO: Store it or write it to a file.

        # Prepare multi-conformer fit.

