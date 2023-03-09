import importlib.resources as files
from string import Template
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="create a file called dftb_in.hsd in the current directory fails if one already exist")
    parser.add_argument(
        '-m',
        '--xtb_method',
        type=str,
        default='GFN2-xTB',
        help="xTB method: ['GFN1-xTB', 'GFN2-xTB', 'IPEA1-xTB']")
    parser.add_argument(
        '-k',
        '--k_pts',
        type=int,
        nargs=3
        default=[2,2,2],
        help="Monkhorst-Pack k-points tuple as list of int will alos automatically compute whether an offset is needed")
    parser.add_argument(
        '-n',
        '--socket_name',
        type=str,
        default="dftbplus",
        help=("name of the unix socket to be used: 'QR', 'DivideAndConquer', 'RelativelyRobust', 'MAGMA'"))
    parser.add_argument(
        '-s',
        '--hamiltonian_solver',
        type=str,
        default=-1.0,
        help="name of the hamiltonian solver")
    parser.add_argument(
        '-t', 
        '--scc_tol',
        type=str, 
        default=1e-7, help="SCC cycle tolerence")

    parser.add_argument(
        '-m', 
        '--max_scc_cycles',
        type=int, 
        default=300, help="Max allowed SCC cycle")
    
    args = parser.parse_args()
    kpts_tuple = [str(kpt) for kpt in args.k_pts]
    offset_tuple = []
    for kpt in args.k_pts:
        if (kpt % 2) == 0:
            offset_tuple.append("0.5")
        else:
            offset_tuple.append("0.0")
    
    variables = {
        "socket_name": args.socket_name,
        "xtb_method": args.xtb_method,
        "max_scc_cycles": args.max_scc_cycles,
        "kpt_a": kpts_tuple[0],
        "kpt_b": kpts_tuple[1],
        "kpt_c": kpts_tuple[2],
        "offset_a": offset_tuple[0],
        "offset_b": offset_tuple[1],
        "offset_c": offset_tuple[2],
        "hamiltonian_solver": args.hamiltonian_solver,
    }

    dftb_in_hsd_template_file = files('tcpg_phonon_tools.templates.dftbp_socket_calc').joinpath('dftb_in.hsd.template').read_text()

    dftb_in_hsd_template = Template(dftb_in_hsd_template_file)

    dftb_in_hsd_result = dftb_in_hsd_template.substitute(variables)

    with open("dftb_in.hsd", "x") as f:
        f.write(dftb_in_hsd_result)




if __name__ == "__main__":
    main()