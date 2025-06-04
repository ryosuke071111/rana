# Coarse-Grained MD Simulation with Martini

This folder contains a minimal example for running a coarse-grained molecular dynamics (MD)
simulation of Bisphenol A diglycidyl ether (DGEBA) and Hexamethylenediamine (HMDA)
using the Martini force field in Python with [OpenMM](https://openmm.org/).

The provided scripts generate a simple box filled with coarse-grained DGEBA and HMDA
molecules (around 100k beads by default) and run a short MD simulation.

## Requirements

* Python 3.8 or later
* `pip`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the simulation with:

```bash
python run_simulation.py
```

The script will output basic energy information to the console and write the final
coordinates to `output.pdb`.

## Notes

The example uses a very simplified mapping and generic Martini-like parameters to
keep the code self‑contained. For production studies you should replace the mapping
and force field details with validated parameters.
