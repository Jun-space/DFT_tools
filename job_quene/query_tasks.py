import pymongo
client_name='tasks'
hostname='ln3'
client = pymongo.MongoClient("mongodb://%s:27018"%hostname)
db = client[client_name]
col = db['vasptasks']

#query={'run_info.path':'/public/home/gouhuiyang/Documents/vasp/MgSiO3_Ca/135GPa/Cmcm/0.875_568/QHA/0.95/ph/001'}

entries=col.find()
job_ids=[]
for num,entry in enumerate(entries):
    job_ids.append(entry['job_id'])
    if entry['job_id'] != num:
        print(entry['job_id'],num,entry['run_info']['path'])
        new={'$set':{'job_id':num}}
        query={'_id':entry['_id']}
        col.update_one(query,new)


#from collections import Counter
#cnt = Counter(job_ids)
#duplicates_with_count = {k: v for k, v in cnt.items() if v > 1}
#print(duplicates_with_count)


#print(col.find_one(query))

#print(col.count_documents(query))

query={'run_info.status':1}
print(col.find_one(query))
#new={'$set':{'run_info.status':0}}
#job_id=2778
#query={'job_id':job_id}
#col.update_one(query,new)

#query={"run_info.path": {'$regex': '/public/home/gouhuiyang/Documents/vasp/MgSiO3_Ca/135GPa/Cmcm/0.666667_89/QHA'}}
#print(col.count_documents(query))
entries=col.find(query)
for entry in entries:
    print(entry['_id'],entry['job_id'],entry['run_info']['path'],entry['run_info']['status'])
    #print(entry['run_info']['path'])
#    new={'$set':{'run_info.status':1}}
#    job_id=entry['job_id']
#    query={'job_id':job_id}
#    col.update_one(query,new)



#new={'$set':{'run_info.status':0}}
#job_id=2778
#query={'job_id':job_id}
#print(col.find_one(query))
#col.update_one(query,new)
#print(col.find_one(query))

#new={'$set':{'run_info.status':0}}
#job_id=2812
#query={'job_id':job_id}
#col.update_one(query,new)

#new={'$set':{'run_info.status':1}}
#job_id=2822
#query={'job_id':job_id}
#col.update_one(query,new)

#new={'$set':{'run_info.status':1}}
#job_id=2826
#query={'job_id':job_id}
#col.update_one(query,new)
