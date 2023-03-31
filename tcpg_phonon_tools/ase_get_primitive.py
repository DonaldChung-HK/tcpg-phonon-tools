#!/usr/bin/env python
import argparse
from spglib import spglib
from ase import Atoms
from ase.io import read


def get_primirive(
        ase_atoms,
        SYMPREC=1e-05,
        ANGLE_TOL=-1.0
    ):
    """Return a ase.Atoms object of a primitive cell found using spglib

    Args:
        ase_atoms (_type_): _description_
        SYMPREC (_type_, optional): _description_. Defaults to 1e-05.
        ANGLE_TOL (float, optional): _description_. Defaults to -1.0.

    Raises:
        RuntimeError: _description_

    Returns:
        _type_: _description_
    """
    cell, positions, atomic_numbers = spglib.find_primitive(
    ase_atoms, symprec=SYMPREC, angle_tolerance=ANGLE_TOL)
    if positions is None:
        print("This space group doesn't have a more primitive unit cell.")
        raise RuntimeError
    else:
        result = Atoms(
                scaled_positions=positions,
                cell=cell,
                numbers=atomic_numbers,
                pbc=True)
    

    lattice_param = result.cell.cellpar()
    print(lattice_param)
    print("Result:")
    print("--------------")
    print("Unit Cell: a: {}, b: {}, c: {}, alpha: {}, beta: {}, gamma: {}".format(*lattice_param))

    return result
        
def main():
    parser = argparse.ArgumentParser(
        description="Find a primitive unit cell using pyspglib and write to a file")
    parser.add_argument(
        'input_file',
        type=str,
        default='POSCAR',
        help="Path to crystal structure file, recognisable by ASE")
    parser.add_argument(
        '--input_format',
        type=str,
        default=None,
        help="Format for input file (needed if ASE can't guess from filename)")
    parser.add_argument(
        '-t',
        '--threshold',
        type=float,
        default=1e-05,
        help=("Distance threshold in AA for symmetry reduction "
              "(corresponds to spglib 'symprec' keyword)"))
    parser.add_argument(
        '-a',
        '--angle_tolerance',
        type=float,
        default=-1.0,
        help="Angle tolerance for symmetry reduction")
    parser.add_argument(
        '-o', '--output_file', default=None, help="Path/filename for output")
    parser.add_argument(
        '--output_format',
        type=str,
        default=None,
        help="Format for input file (needed if ASE can't guess from filename)")
    args = parser.parse_args()

    if args.input_format is None:
        ase_atoms = read(args.input_file)
    else:
        ase_atoms = read(args.input_file, format=args.input_format)

    result = get_primirive(
        ase_atoms,
        SYMPREC=args.threshold,
        ANGLE_TOL=args.angle_tolerance,
    )
    

    if args.output_file == None:
        print('No output file inputed')
        pass
    else:
        if args.output_format == None:
            result.write(args.output_file)
        else:
            result.write(args.input_file, format=args.input_format)
    

    