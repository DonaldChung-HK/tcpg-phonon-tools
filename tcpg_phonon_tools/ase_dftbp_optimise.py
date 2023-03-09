import sys
import argparse
from subprocess import Popen
from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator
from ase.constraints import UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry
from tcpg_phonon_tools.ase_get_primitive import get_primirive


UNIXSOCKET = 'dftbplus'
DFTBP_PATH = 'OMP_NUM_THREADS=1 mpirun -np 1 dftb+'
GEO_PATH = './geo_start.gen'
INPUT_FORMAT = 'gen'
OPTIMIZER = "BFGS"
FMAX = 1.00E-08
FIX_SYM = True
LATTICE_OPT = True
GET_PRIMITIVE = True

SYMPREC = 1e-05
ANGLE_TOL = -1.0

def ase_dftbp_optimise():
    '''Main driver routine.'''

    system = read(GEO_PATH, format=INPUT_FORMAT)

    if GET_PRIMITIVE:
        system = get_primirive(
            system,
            SYMPREC=SYMPREC,
            ANGLE_TOL=ANGLE_TOL,
        )
        
    write('geo.gen', system, format='gen')
    
    if FIX_SYM:
        system.set_constraint(FixSymmetry(system))

    if LATTICE_OPT:
        system_in = UnitCellFilter(system)
    else:
        system_in = system
    
    opt = BFGS(system_in, trajectory='opt.traj', logfile='opt.log')
    with SocketIOCalculator(log=sys.stdout, unixsocket=UNIXSOCKET) as calc:
        Popen([DFTBP_PATH], shell=True)
        system.set_calculator(calc)
        #print("Initial Energy", system.get_potential_energy())
        opt.run(fmax=1.00E-08)

    forces = system.get_forces()
    energy = system.get_potential_energy()
    write('geo_end.gen', system, format='gen')
    print("Final Force", forces)
    print("Final Energy", energy)

def main():
    """run optimisation using dftb+ with the SOCketIO calculator"""
    parser = argparse.ArgumentParser(
        description="Convert between crystal file formats with ASE")
    parser.add_argument('input_file', type=str,
                        help="Path to crystal structure to be converted")
    parser.add_argument('output_file', type=str,
                        help="Output filename")
    parser.add_argument('-f', '--from_format', default=None,
                        help="""Input file format [If this argument is omitted,
                             ASE will guess]""")
    parser.add_argument('-t', '--to_format', default=None,
                        help="""Output file format [If this argument is
                            omitted, ASE will guess]""")
    args = parser.parse_args()


if __name__ == "__main__":
    main()