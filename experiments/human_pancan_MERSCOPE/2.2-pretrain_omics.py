from addict import Dict
from spider import Spider
# os.environ['SLURM_JOB_NAME'] = 'interactive'

model = Spider(work_dir="./results/spider_6_10", recover=True)
model.fit_sc()