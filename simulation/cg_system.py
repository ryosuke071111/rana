import numpy as np
from openmm import unit, openmm
from openmm.app import Topology, Element

# Martini-like bead parameters (simplified)
BEAD_PARAMS = {
    'P1': (0.47 * unit.nanometer, 5.0 * unit.kilojoule_per_mole),
    'P2': (0.47 * unit.nanometer, 5.0 * unit.kilojoule_per_mole),
}

DGEBA_BEADS = ['P1'] * 8  # simple linear model with 8 beads
HMDA_BEADS = ['P2'] * 4   # simple linear model with 4 beads

MASS = 72.0 * unit.dalton
BOND_LENGTH = 0.47 * unit.nanometer
BOND_K = 1250 * unit.kilojoule_per_mole / unit.nanometer ** 2


def build_system(num_dgeba=10000, num_hmda=5000, box_size=50.0):
    """Create topology, system and coordinates for the CG simulation.

    Parameters
    ----------
    num_dgeba : int
        Number of DGEBA molecules.
    num_hmda : int
        Number of HMDA molecules.
    box_size : float
        Simulation box edge length in nanometers.
    """
    topology = Topology()
    system = openmm.System()
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
    nb.setCutoffDistance(1.2 * unit.nanometer)
    bond = openmm.HarmonicBondForce()

    positions = []
    rng = np.random.default_rng()

    def add_molecule(beads):
        chain = topology.addChain()
        residue = topology.addResidue('MOL', chain)
        particle_indices = []
        base = rng.random(3) * box_size
        for bead in beads:
            atom = topology.addAtom(bead, Element.getByAtomicNumber(6), residue)
            system.addParticle(MASS)
            sigma, eps = BEAD_PARAMS[bead]
            nb.addParticle(0.0, sigma, eps)
            pos = base + rng.normal(scale=0.1, size=3)
            positions.append(pos)
            particle_indices.append(len(positions)-1)
        for i in range(len(particle_indices)-1):
            bond.addBond(particle_indices[i], particle_indices[i+1], BOND_LENGTH, BOND_K)

    for _ in range(num_dgeba):
        add_molecule(DGEBA_BEADS)
    for _ in range(num_hmda):
        add_molecule(HMDA_BEADS)

    system.addForce(nb)
    system.addForce(bond)
    system.setDefaultPeriodicBoxVectors(
        [box_size, 0, 0]*unit.nanometer,
        [0, box_size, 0]*unit.nanometer,
        [0, 0, box_size]*unit.nanometer,
    )

    return topology, system, unit.Quantity(np.array(positions), unit.nanometer)
