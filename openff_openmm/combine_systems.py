
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk.openmm.app import PDBFile
import parmed

import time
import sys

ligand_sdf_file = sys.argv[1]  # same atom org as pdb
ligand_pdb_file = sys.argv[2]  # same atom org as sdf
ligand_out_string = sys.argv[3]
receptor_prmtop_file = sys.argv[4]  # parameterized in tleap
receptor_inpcrd_file = sys.argv[5]  # parameterized in tleap
system_out_string = sys.argv[6]

### LIGAND
# OPENFF - load in sdf_file of the ligand for topology
print('Starting to load ligand files.')
start = time.process_time()
ligand_molecule_sdf = Molecule(ligand_sdf_file)
ligand_topology = ligand_molecule_sdf.to_topology()

# OPENFF - load in the FF to be used
force_field = ForceField('openff_unconstrained-1.0.0.offxml')
# OPENFF - couple topology and FF
ligand_system = force_field.create_openmm_system(ligand_topology)     # takes some time

# OPENMM - load in the pdb file of the ligand for coordinates
ligand_pdbfile = PDBFile(ligand_pdb_file)

# PARMED - couple topology, FF, and coordinates to create a parameterized molecule object
ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)  # ligand_pdbfile.topology should be basically empty

# PARMED - save molecule files
ligand_structure.save(ligand_out_string+'.pdb')
ligand_structure.save(ligand_out_string+'.prmtop')

end = time.process_time()
print('Finished handling ligand. Took %.2f seconds'%(end-start))

### RECEPTOR
# PARMED - load in prmtop and inpcrd files into a parmed structure object
print('Starting to load receptors files.')
start = time.process_time()
receptor_structure = parmed.load_file(receptor_prmtop_file,receptor_inpcrd_file)

# PARMED - combining structure objects
complex_structure = receptor_structure + ligand_structure

# PARMED - save molecule files
complex_structure.save(system_out_string+'.pdb')
complex_structure.save(system_out_string+'.prmtop')

end = time.process_time()
print('Finished creating complex structure. Took %.2f seconds'%(end-start))

