
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk.openmm.app import PDBFile
import parmed

import time
import sys
import glob

ligand_smi_file = sys.argv[1]
#ligand_sdf_file = sys.argv[1]  # same atom org as pdb
ligand_pdb_files = glob(sys.argv[2])  # same atom org as sdf
receptor_prmtop_file = sys.argv[3]  # parameterized in tleap
receptor_inpcrd_file = sys.argv[4]  # parameterized in tleap

ligand_pdb_files.sort()

# OPENFF - load in the FF to be used
force_field = ForceField('openff_unconstrained-1.0.0.offxml')

# CREATE SDF FOR LIGANDS - based on the SMILES string
with open(ligand_smi_file,'r') as smi_file:
    for line in smi_file:   # format is SMILES_string molecule_name
        if line[0] == '#':
            continue
        temp = line.split()
        # RDKIT - create sdf file for each ligand
        ligand = Chem.MolFromSmiles(temp[0])
        ligand = Chem.AddHs(ligand)
        writer = SDWriter(temp[1]+'.sdf')
        # OPENFF/OPEN
        lig_pdbs = [pdb for pdb in ligand_pdb_files if pdb.split('.pdb')[0] == temp[1]]   # grabs instances of docked structure to use as the cartesian coords 
        if len(lig_pdbs) > 0:
            temp_pdb = lig_pdbs[0]
            # OPENFF - load in sdf_file of the ligand for topology
            ligand_molecule_sdf = Molecule(temp[1]+'.sdf')
            ligand_topology = ligand_molecule_sdf.to_topology()
            # OPENFF - couple topology and FF
            ligand_system = force_field.create_openmm_system(ligand_topology)     # takes some time
            # OPENMM - load in the pdb files of the ligand for coordinates
            ligand_pdbfile = PDBFile(temp_pdb)
            # PARMED - couple topology, FF, and coordinates to create a parameterized molecule object
            ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)  # ligand_pdbfile.topology should be basically empty
            # PARMED - save molecule files
            ligand_structure.save(ligand_name+'.prmtop')
        else:
            print('No pdbs found for ligand %s; not able to make prmtop file without a pdb'%(temp[1]))
            break

#### LIGAND
## OPENFF - load in sdf_file of the ligand for topology
#print('Starting to load ligand files.')
#start = time.process_time()
#ligand_molecule_sdf = Molecule(ligand_sdf_file)
#ligand_topology = ligand_molecule_sdf.to_topology()
#
## OPENMM - load in the pdb files of the ligand for coordinates
#ligand_pdbfile = PDBFile(ligand_pdb_files[0])
#
## PARMED - couple topology, FF, and coordinates to create a parameterized molecule object
#ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)  # ligand_pdbfile.topology should be basically empty
#
## PARMED - save molecule files
##ligand_structure.save(ligand_name+'.pdb')
#ligand_structure.save(ligand_name+'.prmtop')
#
#end = time.process_time()
#print('Finished handling ligand structures. Took %.2f seconds'%(end-start))
#
#### RECEPTOR
## PARMED - load in prmtop and inpcrd files into a parmed structure object
#print('Starting to load receptors files.')
#start = time.process_time()
#receptor_structure = parmed.load_file(receptor_prmtop_file,receptor_inpcrd_file)
#
#### COMBINED SYSTEM
#for pdb in ligand_pdb_files:
#    ligand_name = pdb.split('.pdb')[0]
#    ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)  # ligand_pdbfile.topology should be basically empty
#    # PARMED - combining structure objects
#    complex_structure = receptor_structure + ligand_structure
#    # PARMED - save molecule files
#    complex_structure.save(system_out_string+'.pdb')
#    complex_structure.save(system_out_string+'.prmtop')
#
#end = time.process_time()
#print('Finished creating complex structure. Took %.2f seconds'%(end-start))

