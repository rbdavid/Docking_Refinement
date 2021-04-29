
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk.openmm.app import PDBFile
import parmed

import sys

sdf_file = sys.argv[1]  # same atom org as pdb
pdb_file = sys.argv[2]  # same atom org as sdf
out_string = sys.argv[3]

# OPENFF - load in sdf_file of the ligand for topology
ligand_molecule_sdf = Molecule(sdf_file)
ligand_topology = ligand_molecule_sdf.to_topology()

# OPENFF - load in the FF to be used
force_field = ForceField('openff_unconstrained-1.0.0.offxml')
# OPENFF - couple topology and FF
ligand_system = force_field.create_openmm_system(ligand_topology)     # takes some time

# OPENMM - load in the pdb file of the ligand for coordinates
ligand_pdbfile = PDBFile(pdb_file)

# PARMED - couple topology, FF, and coordinates to create a parameterized molecule object
ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,ligand_system,xyz=ligand_pdbfile.positions)

# PARMED - save molecule files
ligand_structure.save(out_string+'.pdb')
ligand_structure.save(out_string+'.prmtop')

