import os
import time
from datetime import datetime
import psutil

def process_memory():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss

def profile(*args, **kwargs):
    def wrapper(func):
        obj = kwargs['obj']
        action = func.__name__
        
        tbefore = time.time()
        membefore = process_memory()
        
        result = eval("{obj}.{action}(*args, **kwargs)")
        
        diffMem = process_memory()
        diffTime = time.time() - tbefore
        
        timeExec = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")
        
        folder = kwargs['folder']
        description = kwargs['description']
        
        with open( f'{folder}/log_execution_details.tsv') as f:
            f.write( f"{action}\t{description}\t{timeExec}\t{diffTime}\t{diffMem}\n" )
 
        return result
    return wrapper
