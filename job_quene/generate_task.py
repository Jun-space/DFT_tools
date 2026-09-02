import pymongo
import os
client_name='tasks'
hostname='ln3'
client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
db = client[client_name]
col = db['vasptasks']
start_job_id = db['vasptasks'].count_documents({})
job_id = start_job_id
path=os.getcwd()
run_info = {'slurm_job_id': None, \
                'status': 1, \
                'path': path, \
                'start_time': None, \
                'end_time': None}
vasptask = {'job_id': job_id, \
                'run_info': run_info}
db['vasptasks'].insert_one(vasptask)
