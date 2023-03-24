import subprocess
from string import Template
from pathlib import Path
import warnings
import tcpg_phonon_tools.templates.CRYSTAL.optimisation as optimisation_templates
from importlib.resources import files
from tcpg_phonon_tools.Slurm.SlurmJobHelper import gen_slurm

def opt_crystal_setup(
    label = "foo",
    input_file = "foo.cif",
    basis_set = "POB-TZVP",
    dft_mode = "PBE-D3",
    energy_tol_neg_exp = 7,
    k_pts = (2,2,2),
    nodes = 2,
    wall_time = "10:00:00",
):
    subprocess.run(["cif2cell", "-p", "crystal09", "-o", f"{label}.d12.temp", input_file])
    with open(f"{label}.d12.temp", "r") as f:
        lines = f.readline()

    start = lines.index("CRYSTAL")
    if lines.count("END") == 1 and lines[-1] == "END":
        end = -2 #for geting the second last line
    else:
        warnings.warn("Multiple END or END is not the last line check the cif2cell out file")
        end = -2 #for geting the second last line
    
    geom = '\n'.join(lines[start:end])

    opt_sub_dict = {
        "job_name": "place_holder",
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
        commands=[
            f"mpirun Pcrystal 2> >(tee {label}_opt.out)",
            'find . -name "diis*" -delete',
            'find . -name "fort*" -not -name "fort.34*" -delete'
        ],
    )
