# Collection of My Tools

## Command-Line Arguments
# Positional Argument

`log_file`: Path to the Gaussian log file.
Mutually Exclusive Group

`-opt`: Run in optimization mode to create a JSON file.
`-exc`: Run in excitation energy mode to append excitation energies to an existing JSON file.
`-excopt`: Run in excitation with optimization mode to append both excitation energies and SCF energy for a specified state.
Optional Arguments

`-o, --output`: Specify the output JSON filename. Default is molecule_data.json.
`-state`: Specify the state number for excitation with optimization mode (`-excopt`). Required when using `-excopt`.
Modes of Operation
Optimization Mode (`-opt`)

Extracts molecular charge, multiplicity, optimized geometry, SCF energy, forces, normal modes, and initializes velocities as empty arrays.

Excitation Energy Mode (`-exc`)

Extracts excitation energies and appends them to an existing JSON file.

Excitation with Optimization Mode (`-excopt`)

Extracts both excitation energies and SCF energy for a specified state, appending them to an existing JSON file.

```bash
python parse_gaussian.py molecule.log -opt -o optimized_molecule.json
```

```bash
python parse_gaussian.py molecule.log -exc -o optimized_molecule.json
```

```bash
python parse_gaussian.py molecule.log -excopt -state 1 -o optimized_molecule.json
```

