from pathlib import Path
from subprocess import run
import re
import warnings
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
from importlib.resources import files
import tcpg_phonon_tools.templates.CASTEP as python_script_template_path
import argparse

def phonopy_setup(
    working_dir,
    opt_in_file_name,
    k_pts = (2,2,2),
    supercell = (2,2,2),
    not_gen_slurm = False,
    label = "foo",
    wall_time = "04:00:00",
    nodes = 2,
    path_to_venv = "~/python_env/AMD/bin/activate",
    CASTEP_command = "mpirun castep.mpi "
):
    if not isinstance(working_dir, Path):
        working_dir = Path(str(working_dir))
    
    supercell = ' '.join([str(x) for x in supercell])

    run(['phonopy', '--castep', f'--dim={supercell}', '-d', '-c', opt_in_file_name])

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
        if file.name.startswith('supercell-'): # don't name anything else with supercell in the working folder
            pur = ''.join(re.findall("[0-9]", file.name))
            if pur not in pur_list:
                pur_list.append(pur)
                pur_path = run_path / pur
                pur_path.mkdir(exist_ok=True)
                file.rename(storage_path / f"{pur}.cell")

    pur_list.sort(key=int) #sorting it by number
    
    #verification step to confirm that the displacements id are contineous
    start = 1
    end = len(pur_list) + 1
    test_range = list(range(start, end))
    for i in range(len(pur_list)):
        if test_range[i] != int(pur_list[i]):
            warnings.warn("the displacement list might not be contineous")
    #copy a script from template to folder for customisatin if needed
    castep_python_script_template_file = files(python_script_template_path).joinpath('CastepPhononRunTemplate').read_text()
    with open("run.py", "x") as f:
        f.write(castep_python_script_template_file)
    #setup a slurm to run recurrsively by loading different supercell and running it with the accompanying files
    if not not_gen_slurm:
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
                "Python/3.10.4-GCCcore-11.3.0",
                "CASTEP/21.1.1-iomkl-2021a"
            ],
            set_ups=[
                f"source {path_to_venv}",
                "CASENUM=`printf %03d $SLURM_ARRAY_TASK_ID`",
                f"export CASTEP_COMMAND='{CASTEP_command}'"
            ],
            commands=[
                f"python run.py -f ./storage/$CASENUM.cell -k {k_pts[0]} {k_pts[1]} {k_pts[2]} -p ./run/$CASENUM -l {label}_$CASENUM"
            ]
        )

def main():
    
    parser = argparse.ArgumentParser(
        description="""
        Setup running for CASTEP phonopy run.
        """)
    parser.add_argument('input_file', type=str,
                        help="Input file for phonopy")
    parser.add_argument('-k', '--k_pts',
                        type=int,
                        nargs=3,
                        default=[2,2,2],
                        help="Monkhorst-Pack k-points tuple as list of int will alos automatically compute whether an offset is needed")
    parser.add_argument('-s', '--supercell',
                        type=int,
                        nargs=3,
                        default=[2,2,2],
                        help="diagonal supercell for phonopy")
    parser.add_argument('-n', '--nodes', 
                        default=2,
                        type=int,
                        help="""number of nodes per job""")
    parser.add_argument('-l', '--label', 
                        default="foo",
                        type=str,
                        help="""label name for jobs""")
    parser.add_argument('-t', '--wall_time', 
                        type=str,
                        default="10:00:00",
                        help="""slurm time string for wall time""")
    parser.add_argument('-p', '--path_to_venv', 
                        type=str,
                        default="~/python_env/AMD/bin/activate",
                        help="""path to venv activate script""")
    parser.add_argument('-c', '--CASTEP_command', 
                        type=str,
                        default="mpirun castep.mpi ",
                        help="""command for launching castep""")
    args = parser.parse_args()
    phonopy_setup(
        working_dir = ".",
        opt_in_file_name = args.input_file,
        k_pts = args.k_pts,
        supercell = args.supercell,
        not_gen_slurm = False,
        label = args.label,
        wall_time = args.wall_time,
        nodes = args.nodes,
        path_to_venv = args.path_to_venv,
        CASTEP_command = args.CASTEP_command
    )
