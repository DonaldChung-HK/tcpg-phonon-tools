import sys
import argparse
from subprocess import Popen
from ase import Atoms
from ase.io import read, write

from ase.calculators.socketio import SocketIOCalculator
from ase.constraints import UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from tcpg_phonon_tools.ase_get_primitive import get_primirive

class OptimiserNotIncludedException(Exception):
    "raised when the optimiser is not included"
    pass


def ase_dftbp_optimise(
        input_file = './geo_start.gen',
        input_format = None,
        unix_socket = 'dftbplus',
        dftbp_command = 'dftb+',
        optimiser = "BFGS",
        f_max = 1.00E-08,
        fix_sym = True,
        lattice_opt = True,
        seed_name = "opt"
    ):
    """run a dftbp optimisation using the socket calculator you must use

    Args:
        input_file (str, optional): the path to input file. Defaults to './geo_start.gen'.
        input_format (str or None, optional): the ase format if None it will use the ase input guesser. Defaults to None.
        unix_socket (str, optional): name of the unix socket incase you want to run multiple. Defaults to 'dftbplus'.
        dftbp_command (str, optional): command to run dftb+ you may want to set env variables OMP_STACKSIZE=8G for memory issue and OMP_NUM_THREADS=8 or 1 if not running in parallel. Defaults to 'dftb+'.
        optimiser (str, optional): see the ase.optimizers for details. Defaults to "BFGS".
        f_max (float, optional): the fmax for it to converge. Defaults to 1.00E-08.
        fix_sym (bool, optional): whether to apply the fix symmetry constrain to the system. Defaults to True.
        lattice_opt (bool, optional): whether to optimize cell and position at the same time. Defaults to True.
        seed_name (str, optional): seed name for output files. Defaults to "opt".

    Raises:
        OptimiserNotIncludedException: optimiser not supported in this function

    Returns:
        ase.Atoms: system after optimisation (optimiser should modifly the atoms in place so the result should be included)
    """
    if input_format != None:
        system = read(input_file, format=input_format)
    else:
        system = read(input_file)
        
    write('geo.gen', system, format='gen')
    
    if fix_sym:
        system.set_constraint(FixSymmetry(system))
    
    traj_name = f"{seed_name}.traj"
    log_name = f"{seed_name}.log"


    with SocketIOCalculator(log=sys.stdout, unixsocket=unix_socket) as calc:
        Popen([dftbp_command], shell=True)
        system.set_calculator(calc)
        if lattice_opt:
            filtered = UnitCellFilter(system)
        else:
            filtered = system

        if optimiser == "BFGS":
            from ase.optimize import BFGS
            opt = BFGS(filtered, trajectory=traj_name, logfile=log_name)
        elif optimiser == "LBFGS":
            from ase.optimize import LBFGS
            opt = LBFGS(filtered, trajectory=traj_name, logfile=log_name)
        elif optimiser == "LBFGSLineSearch":
            from ase.optimize import LBFGSLineSearch
            opt = LBFGSLineSearch(filtered, trajectory=traj_name, logfile=log_name)
        elif optimiser == "GPMin":
            from ase.optimize import GPMin
            opt = GPMin(filtered, trajectory=traj_name, logfile=log_name)
        elif optimiser == "FIRE":
            from ase.optimize import FIRE
            opt = FIRE(filtered, trajectory=traj_name, logfile=log_name)
        elif optimiser == "MDMin":
            from ase.optimize import MDMin
            opt = MDMin(filtered, trajectory=traj_name, logfile=log_name)
        else:
            raise OptimiserNotIncludedException(f"{optimiser} is not a recognised optimiser")
        
        opt.run(fmax=f_max)

    forces = system.get_forces()
    energy = system.get_potential_energy()
    system.write(f'{seed_name}_end.gen', format='gen')
    print("Final Force", forces)
    print("Final Energy", energy)
    return system


def main():
    
    parser = argparse.ArgumentParser(
        description="""run optimisation using dftb+ with the SOCketIO calculator \n
    use together with dftbp-socket-setup\n

    to get primitive cell you should run mctools from ajjackson \n
    
    the end trajectory is written to geo_end.gen with trajectory to opt.traj and log to opt.log, \n
    
    the geo.gen is a temperarory file and can be removed afterwards""")
    parser.add_argument('input_file', type=str,
                        help="Path to crystal structure to optimise")
    parser.add_argument('-f', '--input_format', 
                        type=str,
                        default=None,
                        help="""Input file format for ase, if None will use ase.io.read to guess""")
    parser.add_argument('-o', '--optimiser', 
                        type=str,
                        default="BFGS",
                        help="""optimiser for optimisation see ase.optimiser supported: BFGS, LBFGS, LBFGSLineSearch, GPMin, FIRE, MDMin""")
    parser.add_argument('-d', '--dftbp_command', 
                        default="dftb+",
                        type=str,
                        help="""command to run dftbp you can also specify environment variable here such as OMP_STACKSIZE=8G for memory issue and OMP_NUM_THREADS=8 or 1 if not running in parallel""")
    parser.add_argument('-m', '--f_max', 
                        default=1.00E-08,
                        type=float,
                        help="""criteria for the optimiser""")
    parser.add_argument('-u', '--unix_socket', 
                        default="dftbplus",
                        type=str,
                        help="""unix socket to use should match driver.socket.Filename in dftb_in.hsd""")
    parser.add_argument('-n', '--seed_name', 
                        type=str,
                        default="opt",
                        help="""seed name of output file""")
    parser.add_argument('--fix_sym', 
                        action='store_true',
                        help="""fixing symmetry during optimisation""")
    parser.add_argument('--lattice_opt', 
                        action='store_true',
                        help="""optimise the cell and position together""")
    args = parser.parse_args()
    ase_dftbp_optimise(
        optimiser = args.optimiser,
        unix_socket = args.unix_socket,
        dftbp_command = args.dftbp_command,
        input_file = args.input_file,
        input_format = args.input_format,
        f_max = args.f_max,
        fix_sym = args.fix_sym,
        lattice_opt = args.lattice_opt,
        seed_name = args.seed_name     
    )


if __name__ == "__main__":
    main()