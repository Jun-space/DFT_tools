import pymongo
client_name='tasks'
client = pymongo.MongoClient("mongodb://ln3:27018")
db = client[client_name]
col = db['vasptasks']
total=col.count_documents({})
running=col.count_documents({'run_info.status':2})
waiting=col.count_documents({'run_info.status':1})
compelete=col.count_documents({'run_info.status':0})
suberror=col.count_documents({'run_info.status':-2})
print('total tasks: '+str(total))
print('waiting tasks: '+str(waiting))
print('running tasks: '+str(running))
print('complete tasks: '+str(compelete))
print('submit error tasks: '+str(suberror))
