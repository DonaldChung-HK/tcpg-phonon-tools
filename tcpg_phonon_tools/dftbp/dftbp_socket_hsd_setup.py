from importlib.resources import files
from string import Template
import argparse
import tcpg_phonon_tools.templates.dftbp.dftbp_socket_calc as dftbp_socket_calc_template

def gen_dftb_in(
        k_pts = [2,2,2],
        socket_name = "dftbplus",
        max_steps = 1000,
        xtb_method = "GFN2-xTB",
        max_scc_cycles = 300,
        scc_tol = 1e-7,
        hamiltonian_solver = "RelativelyRobust",
        out_file = "dftb_in.hsd"
    ):
    """generate a dftbp input from template. this will provide the basic file for a parameter file to run the socket calculator

    Args:
        k_pts (list(int), optional): kpoints for mp with offset automatically calculated. Defaults to [2,2,2].
        socket_name (str, optional): name of the socket file should match the socket name when runnint the calculation. Defaults to "dftbplus".
        max_steps (int, optional): maximum step for socket calculatoe. Defaults to 1000.
        xtb_method (str, optional): xtb methods 'GFN1-xTB', 'GFN2-xTB', 'IPEA1-xTB'. Defaults to "GFN2-xTB".
        max_scc_cycles (int, optional): max scc per cycle. Defaults to 300.
        scc_tol (float, optional): tollerence for scc. Defaults to 1e-7.
        hamiltonian_solver (str, optional): name of the solver, some might have memory issue 'QR', 'DivideAndConquer', 'RelativelyRobust', 'MAGMA'. Defaults to "RelativelyRobust".
        out_file (str, optional): name of the out file incase you want something else. Defaults to "dftb_in.hsd".
    Returns:
        string: dftbp_in.hsd file content
    """
    kpts_tuple = [str(kpt) for kpt in k_pts]
    offset_tuple = []
    # See SupercellFolding section in DFTB+ manual to check the result
    for kpt in k_pts:
        if (kpt % 2) == 0:
            offset_tuple.append("0.5")
        else:
            offset_tuple.append("0.0")
    
    variables = {
        "socket_name": socket_name,
        "max_steps": max_steps,
        "xtb_method": xtb_method,
        "max_scc_cycles": max_scc_cycles,
        "scc_tol": scc_tol, 
        "kpt_a": kpts_tuple[0],
        "kpt_b": kpts_tuple[1],
        "kpt_c": kpts_tuple[2],
        "offset_a": offset_tuple[0],
        "offset_b": offset_tuple[1],
        "offset_c": offset_tuple[2],
        "hamiltonian_solver": hamiltonian_solver,
    }

    dftb_in_hsd_template_file = files(dftbp_socket_calc_template).joinpath('dftb_in.hsd.template').read_text()

    dftb_in_hsd_template = Template(dftb_in_hsd_template_file)

    dftb_in_hsd_result = dftb_in_hsd_template.substitute(variables)

    with open(out_file, "x") as f:
        f.write(dftb_in_hsd_result)

    return dftb_in_hsd_result

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
        nargs=3,
        default=[2,2,2],
        help="Monkhorst-Pack k-points tuple as list of int will alos automatically compute whether an offset is needed")
    parser.add_argument(
        '-n',
        '--socket_name',
        type=str,
        default="dftbplus",
        help="name of the unix socket file to be used")
    parser.add_argument(
        '-s',
        '--hamiltonian_solver',
        type=str,
        default='RelativelyRobust',
        help="name of the hamiltonian solver: 'QR', 'DivideAndConquer', 'RelativelyRobust', 'MAGMA'")
    parser.add_argument(
        '-t', 
        '--scc_tol',
        type=str, 
        default=1e-7, help="SCC cycle tolerence")
    parser.add_argument(
        '-p', 
        '--max_steps',
        type=int, 
        default=1000, help="Max allowed SCC cycle")
    parser.add_argument(
        '-x', 
        '--max_scc_cycles',
        type=int, 
        default=300, help="Max allowed SCC cycle")
    
    args = parser.parse_args()
    dftb_in = gen_dftb_in(
        k_pts = args.k_pts,
        socket_name = args.socket_name,
        max_steps = args.max_steps,
        xtb_method = args.xtb_method,
        max_scc_cycles = args.max_scc_cycles,
        scc_tol = args.scc_tol,
        hamiltonian_solver = args.hamiltonian_solver,
        out_file = "dftb_in.hsd"
    )
    

if __name__ == "__main__":
    main()