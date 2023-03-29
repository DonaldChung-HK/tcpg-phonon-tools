from pathlib import Path
from subprocess import run
import re
import warnings
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm
from importlib.resources import files
import tcpg_phonon_tools.templates.CRYSTAL.phonopy as phonon_templates
import argparse
from string import Template
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
    copy_input = True,
):
    #loading template
    crystal_INPUT_template_file = files(phonon_templates).joinpath('INPUT.phonopy.template').read_text()
    crystal_INPUT_template = Template(crystal_INPUT_template_file)
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
                file.rename(pur_path / "fort.34")
                if copy_input:
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



    pur_list.sort(key=int) #sorting it by number   
    #verification step to confirm that the displacements id are contineous
    start = 1
    end = len(pur_list) + 1
    test_range = list(range(start, end))
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
