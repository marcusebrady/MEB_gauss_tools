import re
import json
import argparse
from typing import List, Tuple, Dict, Any
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


def atomic_number_to_symbol(atomic_number: int) -> str:
    symbols = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
        7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
        13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
        19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr',
        25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
        31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
        37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
        43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd',
        49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
    }
    return symbols.get(atomic_number, f"Unknown({atomic_number})")


def symbol_to_atomic_number(symbol: str) -> int:
    symbol = symbol.capitalize()
    symbols = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6,
        'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12,
        'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
        'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
        'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
        'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
        'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54,
    }
    return symbols.get(symbol, 0)


def extract_molecular_charge_multiplicity(file_path: str) -> Tuple[int, int]:
    charge = 0
    multiplicity = 1
    with open(file_path, 'r') as file:
        for line in file:
            if 'Charge =' in line and 'Multiplicity =' in line:
                try:
                    parts = line.split()
                    charge_idx = parts.index('Charge') + 2
                    mult_idx = parts.index('Multiplicity') + 2
                    charge = int(parts[charge_idx])
                    multiplicity = int(parts[mult_idx])
                    logging.info(f"Extracted molecular charge: {charge}, multiplicity: {multiplicity}")
                    break
                except (ValueError, IndexError) as e:
                    logging.warning(f"Failed to parse charge and multiplicity from line: {line.strip()} | Error: {e}")
    return charge, multiplicity


atomic_properties = {
    'H': {
        'atomic_number': 1,
        'nuclear_charge': 1,
        'relative_isotopic_mass': 1.008,
        'isotopic_abundance': 99.9885,
        'nuclear_spin': 0.5,
        'nuclear_mass': 1.007825,
    },
    'He': {
        'atomic_number': 2,
        'nuclear_charge': 2,
        'relative_isotopic_mass': 4.0026,
        'isotopic_abundance': 100.0,
        'nuclear_spin': 0.0,
        'nuclear_mass': 4.00260325415,
    },
    'Li': {
        'atomic_number': 3,
        'nuclear_charge': 3,
        'relative_isotopic_mass': 6.94,
        'isotopic_abundance': 7.59,
        'nuclear_spin': 1.0,
        'nuclear_mass': 7.01600455,
    },
    'Be': {
        'atomic_number': 4,
        'nuclear_charge': 4,
        'relative_isotopic_mass': 9.0122,
        'isotopic_abundance': 100.0,
        'nuclear_spin': 2.0,
        'nuclear_mass': 9.0121831,
    },
    'B': {
        'atomic_number': 5,
        'nuclear_charge': 5,
        'relative_isotopic_mass': 10.81,
        'isotopic_abundance': 19.9,
        'nuclear_spin': 3.0,
        'nuclear_mass': 10.0129370,
    },
    'C': {
        'atomic_number': 6,
        'nuclear_charge': 6,
        'relative_isotopic_mass': 12.011,
        'isotopic_abundance': 98.93,
        'nuclear_spin': 0.0,
        'nuclear_mass': 12.0000000,
    },
    'N': {
        'atomic_number': 7,
        'nuclear_charge': 7,
        'relative_isotopic_mass': 14.007,
        'isotopic_abundance': 99.632,
        'nuclear_spin': 1.0,
        'nuclear_mass': 14.0030740,
    },
    'O': {
        'atomic_number': 8,
        'nuclear_charge': 8,
        'relative_isotopic_mass': 15.999,
        'isotopic_abundance': 99.757,
        'nuclear_spin': 0.0,
        'nuclear_mass': 15.99491461956,
    },
}


def get_atomic_properties(symbol: str) -> Dict[str, Any]:
    props = atomic_properties.get(symbol)
    if not props:
        logging.warning(f"Atomic properties for element '{symbol}' not found. Using default values.")
        return {
            'atomic_number': 0,
            'nuclear_charge': 0,
            'relative_isotopic_mass': 0.0,
            'isotopic_abundance': 0.0,
            'nuclear_spin': 0.0,
            'nuclear_mass': 0.0,
        }
    return props


def extract_xyz(file_path: str) -> List[Tuple[int, str, float, float, float]]:
    last_geom: List[Tuple[int, str, float, float, float]] = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if 'Standard orientation:' in line:
                last_geom = []
                start = i + 5
                while start < len(lines):
                    current_line = lines[start].strip()
                    if re.match(r'^-+', current_line):
                        break
                    parts = current_line.split()
                    if len(parts) == 6:
                        try:
                            center_number = int(parts[0])
                            atomic_number = int(parts[1])
                            x = float(parts[3])
                            y = float(parts[4])
                            z = float(parts[5])
                            symbol = atomic_number_to_symbol(atomic_number)
                            last_geom.append((center_number, symbol, x, y, z))
                        except ValueError:
                            logging.debug(f"Skipping line {start+1} due to parsing error: {current_line}")
                            pass
                    start += 1
    if not last_geom:
        raise ValueError("No 'Standard orientation' sections found or no atoms extracted.")
    return last_geom


def extract_scf_energy(file_path: str) -> float:
    scf_en = 0.0
    with open(file_path, 'r') as file:
        for line in file:
            if "SCF Done:" in line:
                try:
                    scf_en = float(line.split()[4])
                    logging.info(f"Extracted SCF energy: {scf_en} Hartrees.")
                except (IndexError, ValueError):
                    logging.warning(f"Failed to parse SCF energy from line: {line.strip()}")
                    scf_en = 0.0
    return scf_en


def extract_forces(file_path: str) -> List[Tuple[int, str, float, float, float]]:
    forces = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in reversed(range(len(lines))):
            if 'Forces (Hartrees/Bohr)' in lines[i]:
                start = i + 3
                while start < len(lines):
                    current_line = lines[start].strip()
                    if re.match(r'^-+', current_line):
                        break
                    parts = current_line.split()
                    if len(parts) == 5:
                        try:
                            center_number = int(parts[0])
                            atomic_number = int(parts[1])
                            fx = float(parts[2]) * 0.529177249
                            fy = float(parts[3]) * 0.529177249
                            fz = float(parts[4]) * 0.529177249
                            symbol = atomic_number_to_symbol(atomic_number)
                            forces.append((center_number, symbol, fx, fy, fz))
                        except ValueError:
                            logging.debug(f"Skipping force line {start+1} due to parsing error: {current_line}")
                            pass
                    start += 1
                break
    if not forces:
        logging.warning("No forces found in the log file.")
    return forces


def extract_freq_info(file_path: str, num_atoms: int) -> List[Tuple[int, str, float, float, float]]:
    normal_modes = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    num_lines = len(lines)
    i = 0
    while i < num_lines:
        line = lines[i]
        if 'Frequencies --' in line:
            freq_values = re.findall(r"[-+]?\d*\.\d+|\d+", line)
            frequencies = [float(f) for f in freq_values]
            num_freqs = len(frequencies)
            if num_freqs == 0:
                logging.warning(f"No frequencies found on line {i+1}")
                i += 1
                continue
            logging.info(f"Found {num_freqs} frequencies on line {i+1}")
            displacement_start = i + 4
            if displacement_start >= num_lines:
                logging.warning(f"Not enough lines after 'Frequencies --' at line {i+1}")
                break
            atom_line = lines[displacement_start].strip()
            if not re.match(r'^Atom\s+AN', atom_line):
                logging.warning(f"Expected 'Atom AN...' line at line {displacement_start+1}, found: '{atom_line}'")
                i += 1
                continue
            displacement_start += 1
            for atom_idx in range(num_atoms):
                if displacement_start >= num_lines:
                    logging.warning(f"Not enough displacement vector lines starting at line {displacement_start+1}")
                    break
                displacement_line = lines[displacement_start].strip()
                if not re.match(r'^\d+', displacement_line):
                    logging.warning(f"Expected atom displacement line at line {displacement_start+1}, found: '{displacement_line}'")
                    break
                parts = displacement_line.split()
                if len(parts) < 2 + 3 * num_freqs:
                    logging.warning(f"Not enough displacement data on line {displacement_start+1}, expected {2 + 3 * num_freqs}, got {len(parts)}")
                    displacement_start += 1
                    continue
                try:
                    index_number = int(parts[0])
                    atomic_number = int(parts[1])
                    symbol = atomic_number_to_symbol(atomic_number)
                    for freq_idx in range(num_freqs):
                        base = 2 + freq_idx * 3
                        if base + 2 >= len(parts):
                            logging.warning(f"Not enough displacement data for mode {freq_idx+1} on line {displacement_start+1}")
                            continue
                        nmx = float(parts[base])
                        nmy = float(parts[base + 1])
                        nmz = float(parts[base + 2])
                        normal_modes.append((index_number, symbol, nmx, nmy, nmz))
                except ValueError:
                    logging.warning(f"Could not parse displacement vector on line {displacement_start+1}")
                displacement_start += 1
            i = displacement_start
        else:
            i += 1
    if not normal_modes:
        logging.warning("No normal modes extracted.")
    return normal_modes


def extract_excitation_energies(file_path: str) -> List[Dict[str, Any]]:
    excitation_energies = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
    last_excitation_section_idx = -1
    for i in reversed(range(len(lines))):
        if 'Excitation energies and oscillator strengths:' in lines[i]:
            last_excitation_section_idx = i
            break
    if last_excitation_section_idx == -1:
        logging.warning("No 'Excitation energies and oscillator strengths:' section found in the log file.")
        return excitation_energies
    for line in lines[last_excitation_section_idx + 1:]:
        stripped_line = line.strip()
        if stripped_line.startswith('Excited State'):
            parts = stripped_line.split()
            if len(parts) < 5:
                continue
            try:
                state_number = int(parts[2].rstrip(':'))
                energy_str = parts[4]
                if energy_str.endswith('eV'):
                    energy_eV = float(energy_str[:-2])
                else:
                    energy_eV = float(energy_str)
                excitation_energies.append({
                    "state_number": state_number,
                    "energy_eV": energy_eV
                })
            except (ValueError, IndexError) as e:
                logging.warning(f"Skipping line due to parsing error: {stripped_line} | Error: {e}")
                continue
    if not excitation_energies:
        logging.warning("No valid excitation energies extracted from the last excitation section.")
    return excitation_energies


def save_geometry_to_json(data: Dict[str, Any], output_filename: str) -> None:
    json_str = json.dumps(data, indent=4)
    def inline_normal_modes(match):
        normal_modes_content = match.group(1)
        inner_pattern = r'\[\n\s*([-\d.eE+]+),\n\s*([-\d.eE+]+),\n\s*([-\d.eE+]+)\n\s*\]'
        inner_replacement = r'[\1, \2, \3]'
        inline_modes = re.sub(inner_pattern, inner_replacement, normal_modes_content)
        return f'"normal_modes": [\n{inline_modes}\n]'
    normal_modes_pattern = r'"normal_modes":\s*\[\n((?:\s*\[[^\]]+\],?\n)+)\]'
    json_str = re.sub(normal_modes_pattern, inline_normal_modes, json_str)
    with open(output_filename, 'w') as json_file:
        json_file.write(json_str)
    logging.info(f"Data has been successfully saved to '{output_filename}' with inline lists.")


def main():
    parser = argparse.ArgumentParser(description="Extract data from Gaussian log files and manage JSON output.")
    parser.add_argument("log_file", help="Path to the Gaussian log file.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-opt", action="store_true", help="Run in optimization mode to create JSON file.")
    group.add_argument("-exc", action="store_true", help="Run in excitation energy mode to append to JSON file.")
    group.add_argument("-excopt", action="store_true", help="Run in excitation with optimization mode to append to JSON file.")
    parser.add_argument("-o", "--output", default="molecule_data.json", help="Output JSON filename (default: molecule_data.json).")
    parser.add_argument("-state", type=int, help="State number for excitation/optimization (required for -excopt).")
    args = parser.parse_args()
    if args.excopt and args.state is None:
        parser.error("The -excopt mode requires the -state argument.")
    if args.opt:
        try:
            charge, multiplicity = extract_molecular_charge_multiplicity(args.log_file)
            geometry = extract_xyz(args.log_file)
            logging.info(f"Extracted geometry for {len(geometry)} atoms.")
            scf_energy = extract_scf_energy(args.log_file)
            if scf_energy != 0.0:
                logging.info(f"SCF energy: {scf_energy} Hartrees.")
            else:
                logging.warning("SCF energy not found in the log file.")
            forces = extract_forces(args.log_file)
            if forces:
                logging.info(f"Extracted forces for {len(forces)} atoms.")
            else:
                logging.warning("No forces extracted.")
            num_atoms = len(geometry)
            normal_modes = extract_freq_info(args.log_file, num_atoms)
            if normal_modes:
                logging.info(f"Extracted {len(normal_modes)} normal modes.")
            else:
                logging.warning("No normal modes extracted.")
            force_dict = {force[0]: [force[2], force[3], force[4]] for force in forces}
            normal_modes_dict: Dict[int, List[List[float]]] = {atom[0]: [] for atom in geometry}
            for mode in normal_modes:
                center_number = mode[0]
                displacement = [mode[2], mode[3], mode[4]]
                if center_number in normal_modes_dict:
                    normal_modes_dict[center_number].append(displacement)
                else:
                    normal_modes_dict[center_number] = [displacement]
            data = {
                "molecule": {
                    "charge": charge,
                    "multiplicity": multiplicity,
                    "atoms": [
                        {
                            "center_number": atom[0],
                            "element_symbol": atom[1],
                            "xyz_coordinates": [atom[2], atom[3], atom[4]],
                            "nuclear_properties": get_atomic_properties(atom[1]),
                            "energy_gradients": force_dict.get(atom[0], [0.0, 0.0, 0.0]),
                            "velocities": [0.0, 0.0, 0.0],
                            "normal_modes": normal_modes_dict.get(atom[0], []),
                        } for atom in geometry
                    ]
                },
                "energy": scf_energy,
                "properties_and_their_derivatives": {
                    "energy": "energy_gradients"
                }
            }
            save_geometry_to_json(data, args.output)
        except Exception as e:
            logging.error(f"An error occurred during optimization mode: {e}")
    elif args.exc:
        try:
            excitation_energies = extract_excitation_energies(args.log_file)
            logging.info(f"Extracted {len(excitation_energies)} excitation energies.")
            if not excitation_energies:
                logging.warning("No excitation energies to add. Exiting.")
                return
            try:
                with open(args.output, 'r') as json_file:
                    existing_data = json.load(json_file)
                logging.info(f"Loaded existing JSON data from '{args.output}'.")
            except FileNotFoundError:
                logging.error(f"The JSON file '{args.output}' does not exist. Please run in optimization mode first.")
                return
            except json.JSONDecodeError as e:
                logging.error(f"Failed to parse JSON file '{args.output}': {e}")
                return
            electronic_states = {}
            for excitation in excitation_energies:
                state_num = excitation["state_number"]
                energy = excitation["energy_eV"]
                electronic_states[f"excitation{state_num}"] = energy
            existing_data["electronic_states"] = electronic_states
            save_geometry_to_json(existing_data, args.output)
        except Exception as e:
            logging.error(f"An error occurred during excitation energy mode: {e}")
    elif args.excopt:
        try:
            excitation_energies = extract_excitation_energies(args.log_file)
            logging.info(f"Extracted {len(excitation_energies)} excitation energies.")
            if not excitation_energies:
                logging.warning("No excitation energies to add. Exiting.")
                return
            scf_energy = extract_scf_energy(args.log_file)
            if scf_energy != 0.0:
                logging.info(f"SCF energy at state {args.state}: {scf_energy} Hartrees.")
            else:
                logging.warning("SCF energy not found in the log file.")
            try:
                with open(args.output, 'r') as json_file:
                    existing_data = json.load(json_file)
                logging.info(f"Loaded existing JSON data from '{args.output}'.")
            except FileNotFoundError:
                logging.error(f"The JSON file '{args.output}' does not exist. Please run in optimization mode first.")
                return
            except json.JSONDecodeError as e:
                logging.error(f"Failed to parse JSON file '{args.output}': {e}")
                return
            electronic_states = {}
            for excitation in excitation_energies:
                state_num = excitation["state_number"]
                energy = excitation["energy_eV"]
                electronic_states[f"excitation{state_num}"] = energy
            existing_data["electronic_states"] = electronic_states
            existing_data[f"SCF_at_state_{args.state}"] = scf_energy
            save_geometry_to_json(existing_data, args.output)
        except Exception as e:
            logging.error(f"An error occurred during excitation with optimization mode: {e}")


if __name__ == "__main__":
    main()


