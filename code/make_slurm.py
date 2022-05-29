import numpy as np
import pandas as pd

pyfile = 157
parallel_num = 70

files = int(np.ceil(pyfile/parallel_num))

for i in range(files):
	start = i*parallel_num+1
	end = min((i+1)*parallel_num, 157)
	with open("stage2_{}.slm".format(i), "a") as f:
		print("#!/bin/bash -l", file = f)
		print("#SBATCH --nodes=1", file = f)
		print("#SBATCH --time=48:00:00", file = f)
		print("#SBATCH --ntasks-per-node=1", file = f)
		print("#SBATCH --mem=60g", file = f)
		print("#SBATCH --mail-type=FAIL", file = f)
		print("#SBATCH --mail-user=he000176@umn.edu", file = f)
		print("#SBATCH --array={}-{}\n".format(start,end), file = f)
		print("cd ~/deepRIV/UKB/code/hdl/combine_p/", file = f)
		print("module load R/3.3.3", file = f)
		print("conda activate tcn\n", file = f)
		print("python3 hdl_stage2_${SLURM_ARRAY_TASK_ID}.py", file = f)
		print("module purge", file = f)


