from pathlib import Path
from string import Template
from subprocess import run
from importlib.resources import files
import tcpg_phonon_tools.templates.dftbp.dftbp_phonon as phonon_template
import re
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
import warnings
import argparse

def gen_dftbp_phonopy_job(
    opt_in_file_name,
    woring_dir = Path("."),
    supercell = (2,2,2),
    nodes = 1,
    wall_time = "00:30:00",
    label = "pho",
    k_pts = [2,2,2],
    max_steps = 1000,
    xtb_method = "GFN2-xTB",
    max_scc_cycles = 300,
    scc_tol = 1e-7,
    hamiltonian_solver = "RelativelyRobust",
    path_to_venv = "~/python_env/AMD/bin/activate",
    python_thread = 2,
    mpi_thread = 1,
    OMP_NUM_THREADS = 8,
    OMP_STACKSIZE = "8G",
    hsd_template_copy = True
):
    supercell_str = ' '.join([str(x) for x in supercell])
    run(['phonopy', '--dftb+', f'--dim={supercell_str}', '-d', '-c', opt_in_file_name])
    file_list = list(woring_dir.iterdir())

    run_path = woring_dir / "run"
    run_path.mkdir(exist_ok=True)

    result_path = woring_dir / "result"
    result_path.mkdir(exist_ok=True)

    gen_dftb_in_phonopy_template(
        k_pts = k_pts,
        max_steps = max_steps,
        xtb_method = xtb_method,
        max_scc_cycles = max_scc_cycles,
        scc_tol = scc_tol,
        hamiltonian_solver = hamiltonian_solver
    )

    pur_list = []
    for file in file_list:
        if file.name.startswith('geo.genS-'): # don't name anything else with supercell in the working folder
            pur = ''.join(re.findall("[0-9]", file.name))
            if pur not in pur_list:
                pur_list.append(pur)
                pur_path = run_path / pur
                pur_path.mkdir(exist_ok=True)
                file.rename(pur_path / file.name)

                if hsd_template_copy:
                    apply_dftb_in_phonopy_template(
                        pur = pur,
                        run_path = run_path,
                        file_name = file.name,
                    )

    pur_list.sort(key=int)
    start = 1
    end = len(pur_list) + 1
    test_range = list(range(start, end))
    for i in range(len(pur_list)):
        if test_range[i] != int(pur_list[i]):
            warnings.warn("the displacement list might not be contineous")    

    gen_slurm(
        slurm_param_list=[
            "-p scarf",
            f"--job-name {label}",
            f"--nodes={str(nodes)}",
            "--exclusive",
            "-C amd",
            f"--time={wall_time}",
            f"--array={start}-{end}"
        ],
        modules=[
            "AMDmodules",
            "dftbp",
        ],
        set_ups=[
            f"source {path_to_venv}",
            "CASENUM=`printf %03d $SLURM_ARRAY_TASK_ID`",
        ],
        commands=[
            f"OMP_STACKSIZE={OMP_STACKSIZE} OMP_NUM_THREADS={OMP_NUM_THREADS} mpirun -np {mpi_thread} dftb+ > {label}_{pur}.out"
        ]
    )

    gen_run_py_phonopy(
        label=label,
        thread=python_thread,
        OMP_NUM_THREADS = OMP_NUM_THREADS,
        OMP_STACKSIZE = OMP_STACKSIZE,
        pur_list = pur_list,
        opt_in_file_name = opt_in_file_name,
        supercell= supercell_str
    )

    

def apply_dftb_in_phonopy_template(
    template_file = "dftb_in.hsd.template",
    pur = "001",
    file_name = "foo",
    run_path = "./run"
):
    if not isinstance(run_path, Path):
        run_path = Path(run_path)

    with open(template_file, "r") as f:
        dftb_template = Template(f.read())
    
    result = dftb_template.substitute({"struct_in_file_name": file_name})

    with open(run_path / pur / "dftb_in.hsd", "w") as f:
        f.write(result)
    

def gen_run_py_phonopy(
    label = "foo",
    thread = 2,
    dftb_command = "dftb+",
    OMP_NUM_THREADS = 8,
    OMP_STACKSIZE = "8G",
    pur_list = [],
    out_file = "run.py",
    supercell = "2 2 2",
    opt_in_file_name = "foo.gen"
):

    variables = {
        "label": label,
        "thread": thread,
        "dftb_command": dftb_command,
        "OMP_NUM_THREADS": OMP_NUM_THREADS,
        "OMP_STACKSIZE": OMP_STACKSIZE,
        "pur_list": str(pur_list),
        "supercell": supercell,
        "opt_in_file_name": opt_in_file_name,
    }

    run_py_template_file = files(phonon_template).joinpath('run.py.template').read_text()

    run_py_template = Template(run_py_template_file)

    run_py_result = run_py_template.substitute(variables)

    with open(out_file, "x") as f:
        f.write(run_py_result)

    return run_py_result



def gen_dftb_in_phonopy_template(
        k_pts = [2,2,2],
        max_steps = 1000,
        xtb_method = "GFN2-xTB",
        max_scc_cycles = 300,
        scc_tol = 1e-7,
        hamiltonian_solver = "RelativelyRobust",
        out_file = "dftb_in.hsd.template"
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

    dftb_in_hsd_template_file = files(phonon_template).joinpath('dftb_in.hsd.template').read_text()

    dftb_in_hsd_template = Template(dftb_in_hsd_template_file)

    dftb_in_hsd_result = dftb_in_hsd_template.safe_substitute(variables)

    with open(out_file, "x") as f:
        f.write(dftb_in_hsd_result)

    return dftb_in_hsd_result

def main():
    parser = argparse.ArgumentParser(
        description="Create a phonopy job setup for quicker run")
    parser.add_argument(
        "input_file",
        type=str,
        help="structure from optimisation"
    )
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
    parser.add_argument(
        '--supercell',
        type=int,
        nargs=3,
        default=[2,2,2],
        help="supercell dim for phonopy")
    parser.add_argument(
        '--nodes',
        type=int, 
        default=1, help="num of nodes for slurm run")
    parser.add_argument(
        '--wall_time',
        type=str, 
        default="00:30:00", help="wall time for slurm run")    
    parser.add_argument(
        '--label',
        type=str, 
        default="pho", help="seed name for run") 
    parser.add_argument(
        '--path_to_venv',
        type=str, 
        default="~/python_env/AMD/bin/activate", help="python venv activate") 
    parser.add_argument(
        '--python_thread',
        type=int, 
        default=2, help="python multiprocess thread for python runner")
    parser.add_argument(
        '--mpi_thread',
        type=int, 
        default=2, help="mpi thread for dftb+ run")
    parser.add_argument(
        '--OMP_NUM_THREADS',
        type=int, 
        default=8, help="OpenMP thread for dftb+ run")
    parser.add_argument(
        '--OMP_STACKSIZE',
        type=str, 
        default="8G", help="OpenMP Stacksize thread for dftb+ run")

    args = parser.parse_args()
    gen_dftbp_phonopy_job(
        opt_in_file_name = args.input_file,
        k_pts = args.k_pts,
        max_steps = args.max_steps,
        xtb_method = args.xtb_method,
        max_scc_cycles = args.max_scc_cycles,
        scc_tol = args.scc_tol,
        hamiltonian_solver = args.hamiltonian_solver,
        supercell = args.supercell,
        nodes = args.nodes,
        wall_time = args.wall_time,
        label = args.label,
        path_to_venv = args.path_to_venv,
        python_thread = args.python_thread,
        mpi_thread = args.mpi_thread,
        OMP_NUM_THREADS = args.OMP_NUM_THREADS,
        OMP_STACKSIZE = args.OMP_STACKSIZE

    )
        