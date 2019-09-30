#!/usr/bin/env python
import os
import re
import time
import csv

def parse_output(filename):
    lines = open(filename,'r').readlines()
    niters = lines[1].split()[-1]
    error = lines[2].split()[-1]
    nopts = lines[3].split()[-1]
    return niters, error, nopts


#testcases = [3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,75,100,125,150,175,200,250,300,350,400,450,500,1000,1500,2000,2500,5000]
testcases = [100,125,150,175,200,250,300,350,400,450,500,1000,1500,2000,2500,5000]
#testcases = [3,4]

out_file = open('results.tsv','w')
tsv_writer = csv.writer(out_file,delimiter='\t')
tsv_writer.writerow(['n','time_used','n_iterations','error','n_opts'])
for case in testcases:
    print(case)
    param = open('param.'+str(case),'w')
    param.write(str(case)+'\n')
    param.close()
    start_time = time.perf_counter()
    os.system('./benchmark_solver.out < param.'+str(case)+' > result.'+str(case))
    end_time = time.perf_counter()
    delta_time = end_time-start_time
    niters, error, nopts = parse_output('result.'+str(case))
    os.system('rm -f  param.'+str(case))
    os.system('rm -f result.'+str(case))
    print(case,delta_time,niters,error,nopts)
    tsv_writer.writerow([case,delta_time,niters,error,nopts])
out_file.close()
    
    
    
