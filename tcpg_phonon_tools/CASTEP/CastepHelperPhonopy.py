from pathlib import Path
from subprocess import run
import re
import warnings
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
from importlib.resources import files
import tcpg_phonon_tools.templates.CASTEP as python_script_template_path


def phonopy_setup(
    working_dir,
    opt_in_file_name,
    k_pts = (2,2,2),
    supercell = "2 2 2",
    not_gen_slurm = False,
    label = "foo",
    wall_time = "04:00:00",
    nodes = 2,
    path_to_venv = "path/to/venv/activate",
    CASTEP_command = "mpirun castep.mpi "
):
    if not isinstance(working_dir, Path):
        working_dir = Path(str(working_dir))
    
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
                pur_path = result_path / pur
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
            modules=[],
            set_ups=[
                f"source {path_to_venv}",
                "CASENUM=`printf %03d $SLURM_ARRAY_TASK_ID`"
                f"export CASTEP_COMMAND='{CASTEP_command}'"
            ],
            commands=[
                f"python run.py -f ./storage/$CASENUM.cell -k {k_pts[0]} {k_pts[1]} {k_pts[2]} -p ./run/$CASENUM -l {label}_$CASENUM"
            ]
        )

