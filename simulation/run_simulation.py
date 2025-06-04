from sys import stdout
from openmm import app, unit, openmm
from cg_system import build_system

# Build system with ~100k beads
# Adjust numbers to get exactly 100k if needed
TOP, SYS, POS = build_system(num_dgeba=10000, num_hmda=5000, box_size=50.0)

integrator = openmm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond, 0.02*unit.picoseconds)
simulation = app.Simulation(TOP, SYS, integrator)
simulation.context.setPositions(POS)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
                                                  potentialEnergy=True,
                                                  temperature=True))

print('Running MD...')
simulation.step(5000)

positions = simulation.context.getState(getPositions=True).getPositions()
with open('output.pdb', 'w') as f:
    app.PDBFile.writeFile(TOP, positions, f)
print('Done. Final coordinates written to output.pdb')
