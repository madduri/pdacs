import os,sys,re,string
from time import time
from optparse import OptionParser
import numpy as np
from scidbpy import interface

#sdb = interface.SciDBShimInterface('http://localhost:8080')
sdb = interface.SciDBShimInterface('http://128.55.57.93:23100')
   
def filetoarr(pathname):
    dirs = pathname.split('/')
    #print dirs
    arrname = dirs[-1].split(".")
    #print arrname
    return "T" + arrname[0]
 	
def main():
    # define options
    parser = OptionParser()
    parser.add_option("-i", "--array1")
    parser.add_option("-j", "--array2")
    parser.add_option("-o", "--output")
    parser.add_option("-c", "--commands")
	
    #sdb.query("create array T100<pid:int64,hid:int64> [i=0:*,10000,0]") 
    # parse
    options, args = parser.parse_args()	

    try:
        #print options
        # retrieve options
        input1   = options.array1
        input2 = options.array2
        output  = options.output
	commands = options.commands
	
	#array1:att1       
	input1_arr = input1.split(':')
	arr1 = filetoarr(input1_arr[0])
	att1 = input1_arr[1]
	#print(arr1,att1)
	
	#array2:att2
	if (input2 != None):
		input2_arr = input2.split(':')
		arr2 = filetoarr(input2_arr[1])
		att2 = input2_arr[1]
	
	cmd_array = commands.split(':')
	print output 
	

	query_res = sdb.new_array()
        
	#print('%s,%s'%input1,input2)
	if (cmd_array[0] == "-a"):
		qry_str = "store(aggregate(" + arr1 + "," + cmd_array[1] + "(" + att1 + ")),{res})"
		print qry_str
	
	if (cmd_array[0] == "-w"):
		if (cmd_array[1] == "filter"):
			if (att1 == ""):
				qry_str = "store(filter(" + arr1+ "," + cmd_array[2] + "),{res})"
				print qry_str
        		else:
				qry_str = "store(project(filter(" + arr1+ "," + cmd_array[2] + ")," + att1 + "),{res})"
				print qry_str
				
	sdb.query(qry_str,res=query_res)		
	
	qry_str2 = "save({res},'%s')"%output
	#print qry_str2
	sdb.query(qry_str2,res=query_res)
	#np = query_res.tolist()
	print query_res.toarray()
	f = open(output,'r');
	print f.read()
	#print np
	sdb.query("remove({res})",res=query_res)

    except Exception, ex:
        print('Error running query.py\n' + str(ex))
    sys.exit(0)
if __name__ == "__main__":
    main()
