import pymongo
client_name='tasks'
hostname='ln3'
client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
db = client[client_name]
col = db['vasptasks']

#query={"run_info.path": {'$regex': '/public/home/gouhuiyang/Documents/vasp/MgSiO3_Ca/135GPa/Cmcm/0.666667_89/QHA'}}

query={'run_info.status':1}
#print(col.count_documents(query))
entries=col.find(query)
job_ids=[]
for entry in entries:
    #print(entry['job_id'],entry['run_info']['path'],entry['run_info']['status'])
    job_ids.append(entry['job_id'])

#for i in range(1921,2221):
for job_id in job_ids:
#    new={'$set':{'run_info.status':1}}
#    job_id=i
    query={'job_id':job_id}
#    col.update_one(query,new)
    col.delete_one(query)
