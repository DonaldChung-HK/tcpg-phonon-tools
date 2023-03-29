import subprocess
from string import Template
from pathlib import Path
import warnings
import tcpg_phonon_tools.templates.CRYSTAL.optimisation as optimisation_templates
from importlib.resources import files
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
import argparse

def opt_crystal_setup(
    label = "foo",
    input_file = "foo.cif",
    basis_set = "POB-TZVP",
    dft_mode = "PBE-D3",
    energy_tol_neg_exp = 7,
    k_pts = (2,2,2),
    nodes = 2,
    wall_time = "10:00:00",
    crystal_command = "mpirun Pcrystal",
    path_to_venv = "~/python_env/AMD/bin/activate",
):
    subprocess.run(["cif2cell", "-p", "crystal09", "-o", f"{label}.d12.temp", input_file])
    with open(f"{label}.d12.temp", "r") as f:
        lines = f.readlines()
    start = lines.index("CRYSTAL\n")
    if lines.count("END\n") == 1 and lines[-1] == "END\n":
        end = -1 #for geting the second last line
    else:
        warnings.warn("Multiple END or END is not the last line check the cif2cell out file")
        end = -1 #for geting the second last line
    
    geom = ''.join(lines[start:end]).rstrip('\n')

    opt_sub_dict = {
        "job_name": f"{label}_opt",
        "geom": geom,
        "basis_set": basis_set,
        "dft_mode": dft_mode,
        "kpoint_string": ' '.join(list(str(x) for x in k_pts)),
        "energy_tol_neg_exp": energy_tol_neg_exp,
    }
    optimisation_template_file = files(optimisation_templates).joinpath('INPUT.optimisation.template').read_text()
    optimisation_template = Template(optimisation_template_file)
    with open("INPUT", "x") as f:
        INPUT_result = optimisation_template.substitute(opt_sub_dict)
        f.write(INPUT_result)
    
    gen_slurm(
        slurm_param_list=[
            "-p scarf",
            f"--job-name {label}_opt",
            f"--nodes={str(nodes)}",
            "--exclusive",
            "-C amd",
            f"--time={wall_time}",
        ],
        modules=[
            "AMDmodules",
            "crystal23",
        ],
        set_ups=[
            f"source {path_to_venv}",
        ],
        commands=[
            f"{crystal_command} 2> >(tee {label}_opt.out)",
            'find . -name "diis*" -delete',
            'find . -name "fort*" -not -name "fort.34*" -delete',
            "#uncomment below for setting up phonopy",
        ],
    )


def opt_setup_cli():
    parser = argparse.ArgumentParser(
        description="""Setup run for CRYSTAL23""")
    parser.add_argument('input_file', type=str,
                        help="Path to crystal structure to optimise")
    parser.add_argument('-b', '--basis_set', 
                        type=str,
                        default="POB-TZVP",
                        help="""Crystal Keyword Basis Set""")
    parser.add_argument('-d', '--dft_mode', 
                        type=str,
                        default="PBE-D3",
                        help="""CRYSTAL keyword DFT mode""")
    parser.add_argument('-c', '--crystal_command', 
                        default="mpirun Pcrystal",
                        type=str,
                        help="""command to run Crystal """)
    parser.add_argument('-t', '--energy_tol_neg_exp', 
                        default=7,
                        type=int,
                        help="""criteria for the optimiser as negative exponential integer""")
    parser.add_argument(
        '-k',
        '--k_pts',
        type=int,
        nargs=3,
        default=[2,2,2],
        help="Monkhorst-Pack k-points tuple as list of int")
    parser.add_argument('-n', '--nodes', 
                        type=int,
                        default=2,
                        help="""number of nodes""")
    parser.add_argument('-w', '--wall_time', 
                        type=str,
                        default="10:00:00",
                        help="""wall time of string""")
    parser.add_argument('-l', '--label', 
                        type=str,
                        default="foo",
                        help="""wall time of string""")
    args = parser.parse_args()

    opt_crystal_setup(
        label = args.label,
        input_file = args.input_file,
        basis_set = args.basis_set,
        dft_mode = args.dft_mode,
        energy_tol_neg_exp = args.energy_tol_neg_exp,
        k_pts = args.k_pts,
        nodes = args.nodes,
        wall_time = args.wall_time,
        crystal_command = args.crystal_command       
    )