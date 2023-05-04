from ase.optimize import BFGS
from ase import units, Atoms
from ase.io import read
import torchani
from ase.spacegroup.symmetrize import FixSymmetry
from ase.constraints import UnitCellFilter
from pathlib import Path
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import numpy as np
import argparse

def torchANI_phonon(input_file,
                    dim = [2, 2, 2],
                    ):
    """This requires the edited phonopy.save function by passing the physical units as well

    Args:
        input (_type_): _description_
        dim (list, optional): _description_. Defaults to [2, 2, 2].
    """
    a = read(input_file)
    a.set_constraint(FixSymmetry(a))
    calculator = torchani.models.ANI2x().double().ase()
    a.calc = calculator
    filtered = UnitCellFilter(a)
    opt = BFGS(filtered, logfile="opt.log", trajectory="opt.traj")
    opt.run(fmax=1e-6)
    a.write("opt_end.gen")


    dim = dim
    atoms_phonopy = PhonopyAtoms(symbols=a.get_chemical_symbols(),
                                scaled_positions=a.get_scaled_positions(),
                                cell=a.get_cell())


    phonopy = Phonopy(atoms_phonopy, supercell_matrix=np.diag(dim),
                    primitive_matrix=None)

    phonopy.generate_displacements()
    supercells = phonopy.supercells_with_displacements
    i = 0
    forces_holder = []
    displacement_dir = Path("disp")
    displacement_dir.mkdir(exist_ok=True)
    for supercell in supercells:

        b = Atoms(
            symbols = supercell.symbols,
            positions = supercell.positions,
            cell = supercell.cell,
        )
        b.calc = calculator
        forces_holder.append(b.get_forces())
        i+=1
        b.write(displacement_dir / f"{i}.gen")
        
    phonopy.forces = forces_holder
    phonopy.produce_force_constants()
    try:
        phonopy.save("phonopy.yaml",settings={'force_constants': True}, physical_units={
                                                                        "length_unit": "angstrom",
                                                                        "force_constants_unit": "eV/angstrom^2"})
    except:
        print("no edited version of phonopy")
        phonopy.save("phonopy.yaml",settings={'force_constants': True})

def torchANI_runner_cli():
    parser = argparse.ArgumentParser(
        description="phono calculation using torchANI")
    parser.add_argument(
        "input_file",
        type=str,
        help="primitive cell in ASE readable format"
    )
    parser.add_argument(
        '-d',
        '--dim',
        type=int,
        nargs=3,
        default=[2,2,2],
        help="supercell dimension diagonal")
    args = parser.parse_args()
    torchANI_phonon(args.input_file,
                    dim = args.dim,
                    )