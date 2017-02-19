#!/usr/bin/env python
#Processes uploads from the user.

# WARNING: Changes in this tool (particularly as related to parsing) may need
# to be reflected in galaxy.web.controllers.tool_runner and galaxy.tools

import urllib, sys, os, gzip, tempfile, shutil, re, gzip, zipfile, codecs, binascii
from scidbpy import interface
from optparse import OptionParser
#from galaxy import eggs
# need to import model before sniff to resolve a circular import dependency
#import galaxy.model
#from galaxy.datatypes.checkers import *
#from galaxy.datatypes import sniff
#from galaxy.datatypes.binary import *

def add_db(filename):
    #sdb = interface.SciDBShimInterface('http://localhost:8080') 
    
    sdb = interface.SciDBShimInterface('http://128.55.57.93:23100') 
    paths = filename.split("/")
    datname = paths[-1].split(".")
    sdb.query("create array T{halo1}<pid:int64,hid:int64> [i=0:*,2000000,0]",halo1=datname[0])
    #print("load(T{halo1},'{filename}')",halo1=datname[0],filename=filename)
    #sdb.query("load(T{halo1},'{filename}')",halo1=datname[0],filename=filename);
    #sdb.query("load(T{halo1},'/home/tmalik/PDACS/galaxy-dist/database/files/000/dataset_121.dat')",halo1=datname[0]);
    str_cmd = 'chmod 777 %s'%filename
    print str_cmd
    os.system(str_cmd) 
    sdb.query("load(T{halo1},'{filename}')",halo1=datname[0],filename=filename);
    return 0

def main():


    parser = OptionParser()
    parser.add_option("-i", "--array1")
    options, args = parser.parse_args()	
    try:
        #print options
        # retrieve options
        input1   = options.array1
        add_db(input1) 
        print "Uploaded succesfully!"
    except Exception, ex:
        print('Error running uploadDB.py\n' + str(ex))
    sys.exit(0)
 
if __name__ == '__main__':
    main()
