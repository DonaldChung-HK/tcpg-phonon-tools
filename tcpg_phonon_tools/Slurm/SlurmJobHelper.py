from pathlib import Path

def gen_slurm(
    output = Path("run.slurm"),
    hashbang = "#! /bin/bash -l",
    slurm_param_list = [
        "-p scarf",
        "--job-name placeholder",
        "--nodes=2",
        "--exclusive",
        "-C amd",
        "--time=00:05:00"
    ],
    modules = [
        "AMDmodules",
        "Python/3.9.6-GCCcore-11.2.0"
    ],
    set_ups = [
      "export FOO=a"  
    ],
    commands = [
        "echo Hello World $FOO"
    ],
):
    if not isinstance(output, Path):
        output = Path(output)

    with open(output, "x") as f:
        #hashbang
        f.write(hashbang + "\n")

        for slurm_param in slurm_param_list:
            param_line = f"#SBATCH {slurm_param}"
            f.write(param_line + "\n")

        #modules lines
        f.write("module purge\n")

        for module in modules:
            module_line = f"module load {module}\n"
            f.write(module_line)

        f.write("module list\n")

        #non_module environment setup
        for set_up in set_ups:
            f.write(set_up + "\n")

        #commands
        for command in commands:
            f.write(command + "\n")


