from addict import Dict
from decipher import DECIPHER
# os.environ['SLURM_JOB_NAME'] = 'interactive'

config = Dict()
config.omics.model.max_steps = 9999999
config.omics.model.epochs = 10
config.device_num = 4


model = DECIPHER(work_dir="./results/decipher_6_10", recover=True)
model.update_config(config)

model.fit_spatial()