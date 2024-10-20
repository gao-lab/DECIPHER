from addict import Dict
from spider import Spider
# os.environ['SLURM_JOB_NAME'] = 'interactive'

config = Dict()
config.omics.model.max_steps = 9999999
config.omics.model.epochs = 10
config.device_num = 4


model = Spider(work_dir="./results/spider_6_10", recover=True)
model.update_config(config)

model.fit_spatial()