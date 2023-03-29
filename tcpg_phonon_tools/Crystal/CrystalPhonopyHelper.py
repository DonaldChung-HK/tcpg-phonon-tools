from pathlib import Path
from subprocess import run
import re
import warnings
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
from importlib.resources import files
import tcpg_phonon_tools.templates.CRYSTAL.phonopy as phonon_templates
import argparse
from string import Template
import shutil
def phonopy_crystal_setup(
    working_dir,
    opt_in_file_name,
    k_pts = (2,2,2),
    supercell = (2,2,2),
    not_gen_slurm = False,
    label = "foo",
    basis_set = "POB-TZVP",
    dft_mode = "PBE-D3",
    energy_tol_neg_exp = 7,
    nodes = 2,
    wall_time = "10:00:00",
    crystal_command = "mpirun Pcrystal",
    path_to_venv = "~/python_env/AMD/bin/activate",
    not_copy_input = False,
    use_symmetry = True,
):
    #loading template
    crystal_INPUT_template_file = files(phonon_templates).joinpath('INPUT.phonopy.template').read_text()
    crystal_INPUT_template = Template(crystal_INPUT_template_file)
    if use_symmetry:
        with open("CRY_SYM", "w") as f:
            f.write("")
    if not isinstance(working_dir, Path):
        working_dir = Path(str(working_dir))
    
    with open(working_dir / "INPUT.template", "w") as f:
        f.write(crystal_INPUT_template_file)

    supercell = ' '.join([str(x) for x in supercell])
    run(['phonopy', '--crystal', f'--dim={supercell}', '-d', '-c', opt_in_file_name])
    file_list = list(working_dir.iterdir())
    run_path = working_dir / "run"
    run_path.mkdir(exist_ok=True)
    #castep gened file are stored as they are only loaded when ase is run
    storage_path = working_dir / "storage"
    storage_path.mkdir(exist_ok=True)

    result_path = working_dir / "result"
    result_path.mkdir(exist_ok=True)
    pur_list = []
    for file in file_list:
        if file.name.startswith('supercell-'):
            if file.name.endswith('.ext'): # don't name anything else with .ext in the working folder
                pur = ''.join(re.findall("[0-9]", file.name))
                if pur not in pur_list:
                    pur_list.append(pur)
                    pur_path = run_path / pur
                    pur_path.mkdir(exist_ok=True)
                    shutil.copyfile(file, storage_path / file.name)
                    file.rename(pur_path / "fort.34")
                    if not not_copy_input:
                        sub_dict = {
                            "job_name": f"{label}_{pur}",
                            "basis_set": basis_set,
                            "dft_mode": dft_mode,
                            "kpoint_string": ' '.join(list(str(x) for x in k_pts)),
                            "energy_tol_neg_exp": energy_tol_neg_exp,
                        }
                        current_INPUT = crystal_INPUT_template.substitute(sub_dict)
                        with open(pur_path / "INPUT", "w") as f:
                            f.write(current_INPUT)
            elif file.name.endswith('.d12'): # this is not usefull but will keep for record
                file.rename(storage_path / file.name)



    pur_list.sort(key=int) #sorting it by number   
    #verification step to confirm that the displacements id are contineous
    start = 1
    end = len(pur_list)
    test_range = list(range(start, end + 1))
    for i in range(len(pur_list)):
        if test_range[i] != int(pur_list[i]):
            warnings.warn("the displacement list might not be contineous")
    
    if not not_gen_slurm:
        gen_slurm(
            slurm_param_list=[
                "-p scarf",
                f"--job-name {label}_pho",
                f"--nodes={str(nodes)}",
                "--exclusive",
                "-C amd",
                f"--time={wall_time}",
                f"--array={start}-{end}"
            ],
            modules=[
                "AMDmodules",
                "crystal23",
                "Python/3.10.4-GCCcore-11.3.0",
            ],
            set_ups=[
                f"source {path_to_venv}",
                "CASENUM=`printf %03d $SLURM_ARRAY_TASK_ID`",
            ],
            commands=[
                "cd ./run/$CASENUM",
                "cp fort.34 fort.34.bak", #crystal will overrite fort.34
                f"{crystal_command} 2> >(tee {label}_$CASENUM.out)",
                f"cp {label}_$CASENUM.out ../../result/",
                'find . -name "diis*" -delete',
                'find . -name "fort*" -not -name "fort.34*" -delete'
            ]
        )

def phonopy_setup_cli():
    parser = argparse.ArgumentParser(
        description="""Setup phonopy for CRYSTAL23""")
    parser.add_argument('input_file', type=str,
                        help="Path to crystal structure optimised")
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
                        default=9,
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
                        help="""wall time  string""")
    parser.add_argument('-l', '--label', 
                        type=str,
                        default="foo",
                        help="""wall time of string""")
    parser.add_argument('--not_copy_input', 
                        action="store_true",
                        help="""not copy the INPUT file template to add arguments""")
    parser.add_argument('--use_symmetry', 
                        action="store_true",
                        help="""use symmetry beware if you have a CRYS_SYM file in the working dir it will use sym anyway""")
    parser.add_argument('--path_to_venv', 
                        type=str,
                        default="~/python_env/AMD/bin/activate",
                        help="""path to venv activate script""")
    parser.add_argument(
        '--supercell',
        type=int,
        nargs=3,
        default=[2,2,2],
        help="diagonal supercell tuple")
    args = parser.parse_args()
    

    phonopy_crystal_setup(
        label = args.label,
        basis_set = args.basis_set,
        dft_mode = args.dft_mode,
        energy_tol_neg_exp = args.energy_tol_neg_exp,
        k_pts = args.k_pts,
        nodes = args.nodes,
        wall_time = args.wall_time,
        crystal_command = args.crystal_command,       
        working_dir = Path("."),
        opt_in_file_name = args.input_file,
        supercell = args.supercell,
        path_to_venv = args.path_to_venv,
        not_copy_input = args.not_copy_input,
        use_symmetry=args.use_symmetry,
    )