import argparse
from pathlib import Path
from typing import Optional
from .input import functionals, methods, geometry, utils

def create_input(args: argparse.Namespace) -> Path:
    """Create ORCA input file from arguments."""
    # Convert SMILES to XYZ if provided
    if args.smiles:
        xyz_path = utils.smiles_to_xyz(args.smiles, Path(args.output_dir) / "molecule.xyz")
    else:
        xyz_path = Path(args.xyz_file)

    # Set up geometry
    geom_builder = geometry.GeometryBuilder(xyz_path)
    geom_lines = geom_builder.create_input_geometry(
        charge=args.charge,
        multiplicity=args.multiplicity
    )

    # Get calculation keywords
    calc_type = methods.CalculationType(args.calc_type.upper())
    functional = functionals.Functional(args.functional.upper())
    
    method_keywords = methods.MethodManager.get_calculation_keywords(calc_type)
    functional_keywords = functionals.FunctionalManager.get_functional_keywords(functional)

    # Combine all input sections
    input_lines = [
        "# ORCA Input File generated by Ophysia",
        *method_keywords,
        *functional_keywords,
        f"! {args.basis}",
        "",
        *geom_lines
    ]

    # Write input file
    output_path = Path(args.output_dir) / f"{xyz_path.stem}.inp"
    with open(output_path, 'w') as f:
        f.write('\n'.join(input_lines))

    return output_path

def main() -> None:
    parser = argparse.ArgumentParser(
        description='OPhysia: ORCA calculation management tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    # Create input file command
    create_parser = subparsers.add_parser('create', help='Create ORCA input file')
    input_group = create_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--smiles', type=str, help='SMILES string of molecule')
    input_group.add_argument('--xyz-file', type=str, help='Path to XYZ file')
    
    create_parser.add_argument('--output-dir', type=str, default='.',
                             help='Output directory')
    create_parser.add_argument('--functional', type=str, default='B3LYP',
                             help='DFT functional')
    create_parser.add_argument('--basis', type=str, default='def2-SVP',
                             help='Basis set')
    create_parser.add_argument('--calc-type', type=str, default='SP',
                             choices=['SP', 'OPT', 'FREQ'],
                             help='Calculation type')
    create_parser.add_argument('--charge', type=int, default=0,
                             help='Molecular charge')
    create_parser.add_argument('--multiplicity', type=int, default=1,
                             help='Spin multiplicity')
    create_parser.add_argument('--n-procs', type=int, default=8,
                             help='Number of processors (requires OpenMPI!)')
    create_parser.add_argument('--submit', action='store_true',
                             help='Submit job after creation')

    args = parser.parse_args()

    if args.command == 'create':
        input_file = create_input(args)
        print(f"Created input file: {input_file} in {args.output_dir}")

if __name__ == "__main__":
    main()