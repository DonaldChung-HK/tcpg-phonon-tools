from datetime import datetime
import ase
import ase.calculators.castep
import ase.io.castep
from ase.io import read
import warnings
from ase.calculators.singlepoint import SinglePointCalculator


def is_fail_in_waring(_warnings):
    for warning in _warnings:
        if "failed" in warning:
            warnings.warn("fail found in Atoms.calc._warnings indicating that geometry might fail to converge")
            return True
    return False


def run_castep_single_point_phonopy(
    system,
    k_pts = (2,2,2),
    on_gamma = True,
    directory = "CASTEP",
    label = "system",
    param_keyword = {
        "task":'SinglePoint',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-7,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
    },
    cell_keyword = {},
    ):
    """wrapper for single point usually using phonopy
    use 'export CASTEP_COMMAND = "mpirun castep.mpi "' tp set the castep command
    param_keyword = {
        "task":'SinglePoint',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-7,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
    }
    cell_keyword = {}
    Args:
        system (ase.Atoms): Atom system to use
        k_pts (tuple, optional): MP kpoints . Defaults to (2,2,2).
        on_gamma (bool, optional): is doing kpoint on gamma. Defaults to True.
        directory (str, optional): name of the directory to put castep in. Defaults to "CASTEP".
        label (str, optional): seed name of the run. Defaults to "system".
        param_keyword (dict, optional): dict of castep keyword for .param .
        cell_keyword (dict, optional): dict of ovcastep keyword for .cell note that the kpoint is set already using calc.set_kpoints.
    """
    # note to self for other keyword
    # geom_modulus_est = 2000 to cheat the optimizer into smalller steps
    # cutoff_energy = int to have a fixed cutoff energy for basis set
    # sedc_apply = False to not apply sedc PBESOL is not compatible with SEDC since the cell is too small

    result = run_castep(
        system = system,
        k_pts = k_pts,
        on_gamma = on_gamma,
        directory = directory,
        label = label,
        param_keyword = param_keyword,
        cell_keyword = cell_keyword,
    )
    return result

def run_castep_opt(
    system,
    k_pts = (2,2,2),
    on_gamma = True,
    directory = "CASTEP",
    label = "system",
    param_keyword = {
        "task":'GeometryOptimization',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-9,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
    },
    cell_keyword = {},
    ):
    """wrapper for geometry optimisation just a wrapper for CASTEP.castep_helper.run_castep() to remember the defaults for geometry optimisation
    use 'export CASTEP_COMMAND = "mpirun castep.mpi "' tp set the castep command
    param_keyword = {
        "task":'GeometryOptimization',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-7,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
    }
    cell_keyword = {}
    Args:
        system (ase.Atoms): Atom system to use
        k_pts (tuple, optional): MP kpoints . Defaults to (2,2,2).
        on_gamma (bool, optional): is doing kpoint on gamma. Defaults to True.
        directory (str, optional): name of the directory to put castep in. Defaults to "CASTEP".
        label (str, optional): seed name of the run. Defaults to "system".
        param_keyword (dict, optional): dict of castep keyword for .param .
        cell_keyword (dict, optional): dict of ovcastep keyword for .cell note that the kpoint is set already using calc.set_kpoints.
    """
    # note to self for other keyword
    # geom_modulus_est = 2000 to cheat the optimizer into smalller steps
    # cutoff_energy = int to have a fixed cutoff energy for basis set
    # sedc_apply = False to not apply sedc PBESOL is not compatible with SEDC since the cell is too small

    result = run_castep(
        system = system,
        k_pts = k_pts,
        on_gamma = on_gamma,
        directory = directory,
        label = label,
        param_keyword = param_keyword,
        cell_keyword = cell_keyword,
    )
    return result

def run_castep_pho(
    system,
    k_pts = (2,2,2),
    on_gamma = True,
    directory = "CASTEP",
    label = "system",
    param_keyword = {
        "task":'PHONON',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-9,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
        "elec_eigenvalue_tol":1e-9,
        "phonon_method":"FINITEDISPLACEMENT",
        "phonon_fine_method":"INTERPOLATE",
    },
    cell_keyword = {
        "phonon_kpoint_mp_grid":"2 2 2",
        "phonon_kpoint_mp_offset":"1/4 1/4 1/4",
    },
    phonon_supercell_matrix = [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
    ]
    ):
    """wrapper for phonon just a wrapper for CASTEP.castep_helper.run_castep() to remember the defaults for phonon and add the supercell matrix to .cell before running
    use 'export CASTEP_COMMAND = "mpirun castep.mpi "' tp set the castep command
    param_keyword = {
        "task":'GeometryOptimization',
        "basis_precision":'precise',
        "xc_functional":'PBE',
        "sedc_scheme":'G06',
        "finite_basis_corr":'AUTOMATIC',
        "geom_force_tol": 1e-4,
        "elec_force_tol": 1e-6,
        "elec_energy_tol": 1e-7,
        "geom_max_iter": 300,
        "grid_scale":2.5,
        "max_scf_cycles":300,
    }
    cell_keyword = {
        "phonon_kpoint_mp_grid":"2 2 2", usually this is 1/2 of the 
        "phonon_kpoint_mp_offset":"1/4 1/4 1/4", formula to set offset is 0.5/phonon_kpoint for on gamma operation
    } 
    Exceptions:
        RuntimeError: raised if dry run fails

    Args:
        system (ase.Atoms): Atom system to use
        k_pts (tuple, optional): MP kpoints . Defaults to (2,2,2).
        on_gamma (bool, optional): is doing kpoint on gamma. Defaults to True.
        directory (str, optional): name of the directory to put castep in. Defaults to "CASTEP".
        label (str, optional): seed name of the run. Defaults to "system".
        param_keyword (dict, optional): dict of castep keyword for .param .
        cell_keyword (dict, optional): dict of ovcastep keyword for .cell note that the kpoint is set already using calc.set_kpoints.
        phonon_supercell_matrix(3x3 list array or None): 3x3 matrix of supercell it will append the cell file since ase castep can't write it if None castep will use 2x2x2
    """
    system = run_castep(
        system = system,
        k_pts = k_pts,
        on_gamma = on_gamma,
        directory = directory,
        label = label,
        param_keyword = param_keyword,
        cell_keyword = cell_keyword,
        not_run=True
    )
    if phonon_supercell_matrix != None:
        with open(f"./{directory}/{label}-out.cell", "a") as f:
            f.write(r"%BLOCK PHONON_SUPERCELL_MATRIX")
            f.write("\n")
            f.write(" {} {} {}\n".format(*phonon_supercell_matrix[0]))
            f.write(" {} {} {}\n".format(*phonon_supercell_matrix[1]))
            f.write(" {} {} {}\n".format(*phonon_supercell_matrix[2]))
            f.write(r"%ENDBLOCK PHONON_SUPERCELL_MATRIX")

    if system.calc.dryrun_ok():
        run_name = str(label)
        potential_energy = system.get_potential_energy()
        phonon_out_cell = read(f"./{directory}/{label}-out.cell")
        holder_calc = SinglePointCalculator(
            phonon_out_cell, 
            energy=potential_energy,
        )
        phonon_out_cell.calc = holder_calc
        print(f"{run_name} is completed at {datetime.now()} with potential energy of : {potential_energy}")
        return phonon_out_cell
    else:
        print("Found error in input")
        print(system.calc._error)
        raise RuntimeError("Found error in input")
        

def run_castep(
    system,
    k_pts = (2,2,2),
    on_gamma = True,
    directory = "CASTEP",
    label = "system",
    param_keyword = {
        "task":'SinglePoint',
    },
    cell_keyword = {},
    not_run = False
    ):
    """generate a normal workflow using ase to run simulation
       use 'export CASTEP_COMMAND = "mpirun castep.mpi "' tp set the castep command
    Args:
        system (ase.Atoms): Atom system to use
        k_pts (tuple, optional): MP kpoints . Defaults to (2,2,2).
        on_gamma (bool, optional): is doing kpoint on gamma. Defaults to True.
        directory (str, optional): name of the directory to put castep in. Defaults to "CASTEP".
        label (str, optional): seed name of the run. Defaults to "system".
        param_keyword (dict, optional): dict of castep keyword for .param .
        cell_keyword (dict, optional): dict of ovcastep keyword for .cell note that the kpoint is set already using calc.set_kpoints.
        not_run(bool, optional): if you need to do some post processing of input file before run this will just setup the run

    Returns:
        ase.Atoms: optimised structure with a dummy calculator attached to it for getting values
    """


    start = datetime.now()

    print(f"Started calculation at {start}")

    print(f"Doing optimisation")

    # to save time copy the castep_keywoard.JSON to the designated place
    # either either one of these path depending on version of ASE
    # ~/.config/ase/castep_keywords.json
    # ~/.ase/castep_keywords.json

    calc = ase.calculators.castep.Castep()
    
    calc._directory = directory

    # include interface settings in .param file
    calc._export_settings = True

    # reuse the same directory
    
    calc._rename_existing_dir = False
    calc._label = label

    # necessary for tasks with changing positions
    # such as GeometryOptimization or MolecularDynamics
    calc._set_atoms = True

    # setting or override any other attribute .param
    for key, value in param_keyword.items():
      setattr(calc.param, key, value)

    # Cell settings
    # calc.cell.kpoint_mp_grid = kpoint_mp_grid
    calc.set_kpts({'size': k_pts, 'gamma': on_gamma})
    
    # setting or override any other attribute in .cell
    for key, value in cell_keyword.items():
      setattr(calc.cell, key, value)
    
    system.calc = calc

    system.calc.initialize()
    if not_run:
        return system
    else:
        # Check for correct input
        if calc.dryrun_ok():
            run_name = str(label)
            potential_energy = system.get_potential_energy()
            # for some reason the atoms are not modified inplace so need to read from cell
            optimised = read(f"{directory}/{label}-out.cell")
            # stop the run if fail to converge
            if is_fail_in_waring(system.calc._warnings):          
                holder_calc = SinglePointCalculator(
                    optimised, 
                    energy=potential_energy
                )
                optimised.calc = holder_calc
                return optimised 
            print(f"{run_name} is completed at {datetime.now()} with potential energy of : {potential_energy}")

            # castep calculator is not great at remembering stuff since there are 2 files to read so it get bugged sometimes
            # beware that if not converged in limit this will cause the simulation to run multiples
            holder_calc = SinglePointCalculator(
                optimised, 
                energy=potential_energy, 
                forces=system.get_forces(), 
                stress=system.get_stress(),
                charges=system.get_charges())
            optimised.calc = holder_calc
            return optimised
        else:
            print("Found error in input")
            print(calc._error)
    
            
