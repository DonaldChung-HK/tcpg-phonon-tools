from pathlib import Path
from subprocess import run
import re
import shutil
import warnings
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm

def phonopy_setup(
    working_dir,
    opt_in_file_name,
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
            pur = re.search("^0*[1-9][0-9]*$", file.name)[-1] #hoping that no changing the filename structure usually the last set of zero paded number is the correct one
            if pur not in pur_list:
                pur_list.append(pur)
                shutil.copyfile(file, storage_path / file.name)

    pur_list.sort(key=int) #sorting it by number
    
    #verification step to confirm that the displacements id are contineous
    start = 1
    end = len(pur_list) + 1
    test_range = list(range(start, end))
    for i in range(len(pur_list)):
        if test_range[i] != int(pur_list[i]):
            warnings.warn("the displacement list might not be contineous")

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
            set_ups=[
                f"soource {path_to_venv}",
                "CASE_NUM=`printf %03d $SLURM_ARRAY_TASK_ID`"
                f"export CASTEP_COMMAND = '{CASTEP_command}'"
            ],
            commands=[
                f"{CASTEP_command} run.py $CASE_NUM"
            ]
        )

