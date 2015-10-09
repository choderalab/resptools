#!/usr/bin/python

#=============================================================================================
# Module for assigning RESP charges using GAMESS.
# Written by John D. Chodera <jchodera@gmail.com> 2007-01-31
#=============================================================================================
# TODO:
# - Look into which point density to use when computing electrostatic potential.  RESP paper uses
# density of 1 point/A^2, but Antechamber seems to use 6 points/A^2.
#=============================================================================================
# Version control information
__version__ = "$Revision: 1.5 $"
# $Source: /share_kubo/Storage/cvsroot/pymol-tools/BatchResp.py,v $
#=============================================================================================
# Imports
import sys
import os
import os.path
import tempfile                    # For managing temporary directories
import openeye.oechem              # Use OpenEye OEChem toolkit
import openeye.oeomega             # Use OpenEye Omega
import exceptions
#=============================================================================================

class BatchResp:
    """
    A class for constructing RESP-derived HF 6-31G* charges using GAMESS and Antechamber.
    
    """

    # Format string used to name conformer directories.
    def __init__(self, gamess_path = None, antechamber_path = None):
        """
        Initialize the RESP module, ensuring that GAMESS and Antechamber can be found.
        """

        # Set variables.
        self.conformer_directory_format_string = "%05d"

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

    def writeGamessInput(self, molecule, directory, calculation = 'optimize'):
        """
        Set up GAMESS calculation to compute electrostatic potential.

        required arguments:
          molecule - the conformer to compute the electrostatic potential for
          directory - the directory to set up the calculatin in

        optional arguments:
          calculation - one of 'optimize' (hours), 'singlepoint' (minutes), or 'test' (seconds) (defaults to 'optimize')

        """

        # GAMESS header blocks.
        # GAMESS is instructed to carry out an RHF claculation with a 6-31G* basis set, after which
        # the electrostatic potential on concentric Connolly surfaces is evaluated.
        gamess_header = { }
        
        # Header for energy minimization.
        gamess_header['optimize'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=OPTIMIZE COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS NPRINT=8 $END
 $SCF    DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $STATPT NSTEP=100 NPRT=-2 NPUN=-2 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
 """

        # Header for no optimization case.
        gamess_header['singlepoint'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=ENERGY COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS NPRINT=8 $END
 $SCF    DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
 """

        # Header for test case, where initial HUCKEL guess is used instead for speed.
        # This should not be used in production.
        gamess_header['test'] = """\
 $CONTRL SCFTYP=RHF EXETYP=RUN RUNTYP=PROP COORD=UNIQUE $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS NPRINT=8 $END
 $SCF    DIRSCF=.TRUE. FDIFF=.FALSE. $END
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

        # Write input file to work directory.
        self.writeFile(os.path.join(directory, 'gamess.inp'), gamess_input)


    def findLine(self, lines, string):
        """\
        Find the first line that begins with the given string.

        required arguments:
          lines - the lines to search
          string - the given string to check for

        returns:
          index - the index of the line (starting from zero) where the match was found, or 'None' if no match found

        """

        index = 0
        for line in lines:
            if line.startswith(string):
                return index
            index += 1

        raise Exception("Line not found.")   
        #return None          

    def readGamessOutput(self, directory, molecule):
        """\
        Process GAMESS output to extract optimized coordinates, Mulliken charges, and electrostatic potential.

        required arguments:
          directory - directory containing GAMESS output as gamess.log and gamess.dat
          molecule - OEMol molecule

        return values:
          optimized_molecule - the optimized conformation, with Mulliken charges
          potential - the electrostatic potential felt by a test positive charge, where potential[(x,y,z)] = psi, with distance in borhs
          
        """

        print "Reading GAMESS output from %s..." % directory
        
        # Read in file contents.
        log_lines = self.readFile(os.path.join(directory,'gamess.log'))
        dat_lines = self.readFile(os.path.join(directory,'gamess.dat'))

        # Determine number of atoms.
        natoms = molecule.NumAtoms()

        # Extract electrostatic potential.
        import re
        index = self.findLine(log_lines, ' NUMBER OF POINTS SELECTED FOR FITTING =')
        npoints = int(re.match(' NUMBER OF POINTS SELECTED FOR FITTING =\s*(\d+)\n', log_lines[index]).groups()[0])
        index = self.findLine(dat_lines, ' ELECTRIC POTENTIAL, IPT,X,Y,Z,ELPOTT')
        potential = { } # Electrostatic potential: potential[(x,y,z)] = psi
        for line in dat_lines[index+1:index+npoints+1]:            
            # Parse line to extract coordinate triple (x,y,z) and potential psi in bohr units.
            # Fortran format: (1X,I5,2X,3F10.5,E15.6)
            # Example format: 
            #     1     1.59920   3.28673   1.32987  -0.990640E-02
            x = float(line[8:18])
            y = float(line[18:28])
            z = float(line[28:38])
            psi = float(line[38:54])
            
            # Store potential at this grid point.
            potential[(x,y,z)] = psi
        # Check number of points read.
        if (len(potential.keys()) != npoints):
            raise Exception( "Only read %d points (expected %d)." % (len(potential.keys()), npoints) )
            
        # Make a deep copy of the molecule.
        optimized_molecule = openeye.oechem.OEMol(molecule)
        
        # Extract coordinates of optimized geometry.
        # Example format:
        # C1          6.0   0.5931344669  -0.3454239446   0.7037368271        
        index = self.findLine(dat_lines, '----- RESULTS FROM ')
        # TODO: Check if geometry search converged: "SUCCESSFUL" vs "INCOMPLETE"
        coordinates = { } # Atomic coordinates: coords[name] = (x, y, z)
        for line in dat_lines[index+5:index+5+natoms]:
            # Parse line to extract atom name (name) and coordinate triple (x,y,z) and convert to bohrs.
            (name, charge, x, y, z) = line.split()
            x = float(x)
            y = float(y)
            z = float(z)
            
            # Store coordinates in a dictionary by atom name.
            coordinates[name] = (x, y, z)

        # Assign coordinates to atoms in mol2 file.
        for atom in optimized_molecule.GetAtoms():
            optimized_molecule.SetCoords(atom, coordinates[atom.GetName().upper()])

        # Extract Mulliken charges.
        # Example format:
        #          1         2         3         4         5
        #012345678901234567890123456789012345678901234567890
        #C1           6.16547  -0.16547   6.10684  -0.10684        
        index = self.findLine(dat_lines, ' POPULATION ANALYSIS')
        mulliken_charges = { } # Mulliken charges: mulliken_charges[name] = charge
        for line in dat_lines[index+1:index+1+natoms]:
            # Parse line.
            name = line[0:10].strip()
            mulliken_population = float(line[10:20])
            mulliken_charge = float(line[20:30])
            lowdin_population = float(line[30:40])
            lowdin_charge = float(line[40:50])
                                                  
            # Store Mulliken charge.
            mulliken_charges[name] = mulliken_charge
        
        # Assign charges to atoms in mol2 file.
        for atom in optimized_molecule.GetAtoms():
            atom.SetPartialCharge(mulliken_charges[atom.GetName().upper()])

        # Extract energy.
        # Example format:
        #          1         2         3         4         5
        #012345678901234567890123456789012345678901234567890123456789
        #ENERGY IS     -798.2033584902 E(NUC) IS     1054.5798999873
        index = self.findLine(dat_lines, 'ENERGY IS ')
        hartree_to_kcal_per_mole = 627.509391
        line = dat_lines[index]
        energy = float(line[10:29]) * hartree_to_kcal_per_mole
        optimized_molecule.SetEnergy(energy)
                
        return (optimized_molecule, potential)
                    
    def readMol2File(self, filename):
        """\
        Read a mol2 file as an OEChem molecule.

        required arguments:
          filename - the name of the mol2 file
         
        return values:
          molecule - the molecule read in
          
        """

        # Open input mol2 file.
        input_molecule_stream = openeye.oechem.oemolistream()

        # Direct oemolistream to read all molecules in file into a multi-conformer molecule.
        input_molecule_stream.SetConfTest( openeye.oechem.OEIsomericConfTest() )

        input_molecule_stream.open(filename)

        # Create a molecule.
        molecule = openeye.oechem.OEMol()

        # Read the molecule from file.
        openeye.oechem.OEReadMolecule(input_molecule_stream, molecule)

        # Return the moelecule.
        return molecule
    
    def writeMol2File(self, molecule, filename):
        """\
        Write an OEChem molecule to a mol2 file.

        required arguments:
          molecule - an OEChem molecule
          filename - the filename to write to.

        """
        
        # Create an output stream.
        output_molecule_stream = openeye.oechem.oemolostream()

        # Associate it with the given filename.
        output_molecule_stream.open(filename)

        # Write the molecule to the output stream.
        openeye.oechem.OEWriteMolecule(output_molecule_stream, molecule)

        # Close the output stream.
        output_molecule_stream.close()
        
    def setupJobs(self, molecule, molecule_name, directory, generate_conformers = False):
        """\
        Set up GAMESS jobs for computing charges on the QB3 cluster in given subdirectory.

        required arguments:
          molecule - an OEChem molecule
          directory - the name of the directory to set up

        optional arguments:
          generate_conformers - if True, multiple conformers will be generated by OpenEye's Omega (default False)

        returns:
          nconformers - the number of conformers


        """

        print "Setting up jobs for molecule %s in directory %s" % (molecule_name, directory)

        # Generate multiple conformers if requested.
        if (generate_conformers):
            # Initialize an instance of Omega.
            omega = openeye.oeomega.OEOmega()
            # Generate multiple conformers from the given molecule.
            omega(molecule)            

        # Write the molecule.
        self.writeMol2File(molecule, os.path.join(directory, '%s.mol2' % molecule_name))

        # Perform electronic structure calculations on each conformer.
        nconformers = 0 # Count of current number of conformers.
        for conformer in molecule.GetConfs():
            # Increment count of the current number of conformers.
            nconformers += 1

            # Form conformer directory name.
            conformer_directory = os.path.join(directory, self.conformer_directory_format_string % nconformers)

            # Create a directory for this conformation.
            os.mkdir(conformer_directory)

            # Write the conformer .mol2 file to this directory.
            self.writeMol2File(conformer, os.path.join(conformer_directory, '%s.mol2' % molecule_name))            
            
            # Write GAMESS input files.
            self.writeGamessInput(conformer, conformer_directory, calculation = 'optimize')

        # Create an array job input file.
        max_job_length = 10.0 # ten hours

        # Select queue type.
        if (max_job_length >= 12.0):
            queue_type='long.q'
        elif (max_job_length >= 0.5) :
            queue_type='medium.q'
        else:
            queue_type='short.q'

        # Form job script.
        conformer_directory_format_string = self.conformer_directory_format_string
        sge_job_script = """\
#!/bin/tcsh
#$ -l opt64=true
#$ -p -4
#$ -S /bin/tcsh
#$ -cwd
#$ -N %(molecule_name)s
#$ -r y
#$ -j y
#$ -l netapp=1G,scratch=1G,mem_total=3G
#$ -l h_rt=%(max_job_length)s:00:00
#$ -q %(queue_type)s
#$ -t 1-%(nconformers)d

set echo on

date
hostname
pwd

# Run GAMESS in conformer directory.
echo Computing electrostatic potential for conformer $SGE_TASK_ID
cd `echo $SGE_TASK_ID | awk '{printf("%(conformer_directory_format_string)s",$1)}'`
gms gamess.inp -l gamess.log

# report ending date
date

""" % vars()

        # Write the SGE job script
        sge_job_script_filename = os.path.join(directory, 'run_batch.cmd')
        outfile = open(sge_job_script_filename, 'w')        
        outfile.write(sge_job_script)
        outfile.close()

        # Return number of conformers
        return nconformers

    def extractAtomBlockFromRespInputFile(self, filename):
        """
        Extract the "atom" block from a RESP input file.

        required arguments:
          filename - the name of the RESP input file to read

        return values:
          atom_block - an array of the atom block from the RESP file

        """

        # Read file contents into a list of lines.
        lines = self.readFile(filename)

        # Determine which line 'Resp charges for organic molecule' falls on.
        index = self.findLine(lines[1:], 'Resp charges for organic molecule') + 1

        # Extract number of atoms and net charge.
        natoms = int(lines[index+1][5:10])

        # Select atom block.
        atom_block = lines[index+2:index+2+natoms]

        # Return the atom block.
        return atom_block

    def writeRespInputFile(self, filename, molecule, net_charge, atom_block, stage):
        """
        Write a multi-conformer RESP input file.

        required arguments:
          filename - the name of the RESP input file to write
          molecule - the molecule to write
          net_charge - the (integral) net charge of the molecule
          atom_block - list of atomic charges and equivalence classes in (I5,I5) format          
          stage - set to 1 for stage 1, 2 for stage 2 RESP fitting.
          
        """

        # Open file for writing.
        outfile = open(filename, 'w')

        # Write title.
        outfile.write("RESP fit for %d conformers of %s.\n" % (molecule.NumConfs(), molecule.GetTitle()) )

        # Determine number of conformers.
        nconformers = molecule.NumConfs()

        # Determine number of atoms.
        natoms = molecule.NumAtoms()
        
        # Write header.
        if (stage == 1):
            header = """\
 &cntrl

 nmol = %(nconformers)d,
 iunits = 0,
 ihfree = 1,
 ioutopt = 1,

 &end\n""" % vars()
        elif (stage == 2):
            header = """\
 &cntrl

 nmol = %(nconformers)d,
 iunits = 0,
 ihfree = 1,
 ioutopt = 1,
 iqopt = 2,
 qwt = 0.001,

 &end\n""" % vars()
        else:
            raise ParameterException("Unrecognized stage argument: stage = %d\n" % stage)
        
        outfile.write(header)
        
        # Write conformers and weights.
        conformer_index = 0
        for conformer in molecule.GetConfs():
            # Construct conformer directory name.
            conformer_index += 1
            conformer_name = self.conformer_directory_format_string % conformer_index
            
            # Set weight for this conformer.
            # TODO: Maybe this should be based on the GAMESS energy for this conformation?
            weight = 1.0
            
            # TODO: Append geometry and electrostatic potential to RESP input file.
            outfile.write('%10.5f\n' % weight) # weight for this conformer
            outfile.write('%10s\n' % conformer_name) # name of this conformer (irrelevant)
            outfile.write('%5d%5d\n' % (net_charge, natoms)) # net integral charge and number of charge centers
            for line in atom_block:
                outfile.write(line) # nuclear charge and equivalence class in (I5,I5) format
            if(conformer_index < nconformers):
                outfile.write('\n') # blank line between conformers

        # Write charge constraints.
        outfile.write('\n') # none
        outfile.write('\n')
            
        # Write atom equivalences between conformations in (16I5) format (see AMBER8 manual, pg. 235).
        if (nconformers > 1):
            atom_index = 0
            for atom in molecule.GetAtoms():
                atom_index += 1
                outfile.write('%5d\n' % nconformers)
                pairs_written_this_line = 0
                for conformer_index in range(1,nconformers+1):
                    outfile.write('%5d%5d' % (conformer_index, atom_index))
                    pairs_written_this_line += 1
                    
                    if(pairs_written_this_line == 8):
                        outfile.write('\n')
                        pairs_written_this_line = 0
                outfile.write('\n')
                                
        outfile.write('\n')

        # Close file.
        outfile.close()               
    
    def respFit(self, directory, net_charge = 0):
        """
        Perform RESP fit on the molecule in a given directory after completion of GAMESS jobs.

        required arguments:
          directory - the directory to process

        optional arguments:
          net_charge - integral net charge of the molecule

        """

        # Determine name of mol2 file.
        import glob
        molecule_filename = glob.glob(os.path.join(directory,'*.mol2'))[0]
        print "molecule filename = %s" % molecule_filename

        # Extract the name of the mol2 file.
        molecule_name = os.path.basename(molecule_filename).replace('.mol2','')
        print "molecule name = %s" % molecule_name

        # Load the (multiconformer) mol2 file.        
        molecule = self.readMol2File(molecule_filename)

        # Create an Antechamber file from a single conformer mol2 file.
        conformer_directory = os.path.join(directory, self.conformer_directory_format_string % 1)
        mol2_filename = os.path.join(conformer_directory, '%s.mol2' % molecule_name)
        ac_filename = os.path.join(directory, '%s.ac' % molecule_name)
        os.system('antechamber -i %(mol2_filename)s -fi mol2 -o %(ac_filename)s -fo ac -nc %(net_charge)d' % vars() )

        # Use Antechamber's 'respgen' to generate template single-conformer RESP input files for stage 1 and stage 2 fitting.
        resp1_filename = os.path.join(directory, 'resp1.in')
        os.system('respgen -i %(ac_filename)s -o %(resp1_filename)s -f resp1' % vars() )
        resp1_contents = self.readFile(resp1_filename)
        
        resp2_filename = os.path.join(directory, 'resp2.in')
        os.system('respgen -i %(ac_filename)s -o %(resp2_filename)s -f resp2' % vars() )
        resp2_contents = self.readFile(resp2_filename)
    
        # Construct multi-conformer RESP input files.
        atom_block_1 = self.extractAtomBlockFromRespInputFile(resp1_filename)
        atom_block_2 = self.extractAtomBlockFromRespInputFile(resp2_filename)

        # Create resp stage 1 and 2 files.
        self.writeRespInputFile(resp1_filename, molecule, net_charge, atom_block_1, 1)
        self.writeRespInputFile(resp2_filename, molecule, net_charge, atom_block_2, 2)
        
        # Write potential in RESP format.
        esp_filename = os.path.join(directory, 'esp.dat')
        espfile = open(esp_filename, 'w')

        # Process each conformer directory to extract optimized geometry and electrostatic potential.
        print "Reading optimized geometries and electrostatic potentials from GAMESS output and writing potential file..."
        conformer_index = 0
        for conformer in molecule.GetConfs():
            # Construct conformer directory name.
            conformer_index += 1
            conformer_directory = os.path.join(directory, self.conformer_directory_format_string % conformer_index)
            
            # Read GAMESS output files to obtain final optimized geometry and electrostatic potential.
            (optimized_conformer, potential) = self.readGamessOutput(conformer_directory, conformer)

            # Write to ESP file.
            npoints = len(potential.keys())
            espfile.write('%5d%5d\n' % (optimized_conformer.NumAtoms(), npoints))
            # Write atom coordinates in Fortran (17x,3e16.7) format (in bohrs).
            for atom in optimized_conformer.GetAtoms():
                # Convert from Angstroms to Bohr
                angstroms_per_bohr = 0.52917725            
                point = ()
                for coord in optimized_conformer.GetCoords(atom):
                    point += ( coord / angstroms_per_bohr, )
                # Write
                espfile.write('%17s%16.7e%16.7e%16.7e\n' % (('',) + point) )
            # Write potential in Fortran (qpot,x,y,z) in (1x,4e16.7) format
            for point in potential.keys():
                espfile.write(' %16.7e%16.7e%16.7e%16.7e\n' % ((potential[point],) + point) )
        espfile.close()

        # Perform RESP fit.
        # Command-line options employed here (for future reference):
        # -O : overwrite output files
        # -i [resp input file]
        # -o [resp output file]
        # -p [point charge optimization summary]
        # -t [final point charges]
        # -q [initial point charges]
        # -e [input (potentially multiconformer) electrostatic potential file]
        # -s [computed electrostatic potential from optimized point charges?]
        print "Performing RESP fit..."
        os.system('cd %(directory)s ; resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -e esp.dat -s resp1.esp' % vars())
        os.system('cd %(directory)s ; resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg -t resp2.chg -e esp.dat' % vars())

        # Read RESP charges from final charge file 'resp2.chg'.
        charge_lines = self.readFile(os.path.join(directory,'resp2.chg'))
        charge_list = [ ]
        for line in charge_lines:
            for charge in line.split():
                charge_list.append( float(charge) )

        # Read the RESP charges into molecule.
        charged_molecule = openeye.oechem.OEMol(molecule)

        # Store RESP charges in charged molecule.
        atom_index = 0
        for atom in charged_molecule.GetAtoms():
            atom_index += 1
            atom.SetPartialCharge(charge_list[atom_index-1])

        # Write molecule with updated charges to mol2 file.
        charged_molecule_filename = '%s-resp.mol2' % molecule_name
        self.writeMol2File(charged_molecule, charged_molecule_filename)
    
        # Return the charged molecule
        return charged_molecule
