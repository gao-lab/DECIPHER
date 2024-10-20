from addict import Dict
from decipher import DECIPHER
# os.environ['SLURM_JOB_NAME'] = 'interactive'

model = DECIPHER(work_dir="./results/decipher_6_10", recover=True)
model.fit_sc()