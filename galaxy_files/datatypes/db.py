"""
SQLite3 database data type
"""

import pkg_resources
import logging, sqlite3
from data import Data
from galaxy import util
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from sniff import *
import json
import sys

def is_sqlite_3_file(filename):
    txt = open(filename,"rb").read(15)
    return txt == "SQLite format 3"

class Db( Data ):
    file_ext = "db"

    MetadataElement( name="table_names", default=0, desc="Table names", readonly=True, optional=True, visible=True, no_value=0 )
    MetadataElement( name="config", default=0, desc="Configure Parameters", readonly=True, optional=True, visible=True, no_value=0 )
    MetadataElement( name="column_names", default=0, desc="Names of columns in all tables", readonly=True, optional=True, visible=True, no_value=0)

    def get_master_data(self,con, sm):
        cur = con.cursor()
        try:
            # columns names are: name, type, tbl_name, rootpage, sql
            cur.execute(sm)
            rec = cur.fetchall()
            return rec
        except sqlite3.DatabaseError as e:
            return None

#    def get_visualizations( self, dataset ):
#        vizs = []
#        vizs = super( Db, self ).get_visualizations( dataset )
#        vizs.append( 'scat' )
        #vizs.append( 'google.com' )
#        return  vizs

    def to_json_string(self, dataset ):
        # print "hooray! got a dataset"
        # cols = ['a','b'] or the list of column headings
        # rec = select * from table_we_want
        # json.dumps(dict(zip(cols,zip(*rec))))
        con = sqlite3.connect(dataset.file_name)
        cur = con.cursor()
        metarec = self.get_master_data(con, 'select name,sql from sqlite_master where type = "table" and name!="metadata"')
        #metarec = self.get_master_data(con, 'select name from sqlite_master where type = "table" and name!="metadata"')
        # each list entry will be a tuple of (table_name, create_table_statement)
        column_names={}
        col_names=[]
        #pat=re.compile(r"\(([^\)]+)\)")
        for table_name, sql in metarec:
            # clean out the formatting junk in sql
            #sql = sql.replace('\n','')
            #print sql
	    #cols = pat.findall(sql)[0].split()[0::2]
            #print table_name
            q1 = 'PRAGMA table_info(%s)'
            col_rec = self.get_master_data(con, q1%(table_name))
            #print col_rec
            for d in col_rec:
		col_names.append(d[1])
	    column_names[table_name] = col_names
            col_names=[]
            #print column_names[table_name]
        #print "----\n" 
        query = 'select * from %s'
        all={}
        for k,v in column_names.iteritems():
            rec = self.get_master_data(con, query%(k))
            all[k] = dict(zip(v,zip(*rec)))
        #print all
        return json.dumps(all)
        
    # sniff is an instance function only because Galaxy requires it to
    # be.
    def sniff(self, filename):
        return is_sqlite_3_file(filename)

    def set_meta(self, dataset, **kwd):
        print "----\n"
        print "In Db.set_meta"
        print "self is: %s" % str(self)
        print "class of dataset is: %s" % str(dataset.__class__)
        print "dataset is: %s" % str(dataset.__dict__)
        print "----\n"
        con = sqlite3.connect(dataset.file_name)
        metarec = self.get_master_data(con, 'select name,sql from sqlite_master where type = "table" and name!="metadata"')
        dataset.metadata.column_names = {}
        dataset.metadata.table_names = []
        col_names = []
        #pat=re.compile(r"\(([^\)]+)\)")
        for table_name, sql in metarec:
            dataset.metadata.table_names.append(table_name)          
            #sql = sql.replace('\n','')
            #cols = pat.findall(sql)[0].split()[0::2]
            #dataset.metadata.column_names[table_name] = cols
            #print dataset.metadata.column_names[table_name]
            q1 = 'PRAGMA table_info(%s)'
            col_rec = self.get_master_data(con, q1%(table_name))
            #print col_rec
            for d in col_rec:
                col_names.append(d[1])
            dataset.metadata.column_names[table_name] = col_names
            col_names=[]
            print dataset.metadata.column_names[table_name]
 
        rec = self.get_master_data(con, 'select * from metadata')
        y=""
        for n in rec:
            y+="%s=%s\n"%n
        dataset.metadata.config = y
        return True

    def get_meta(self, dataset, **kwd):
	con = sqlite3.connect(dataset.file_name)
        metarec = self.get_master_data(con, 'select name,sql from sqlite_master where type = "table" and name!="metadata"')
        column_names = {}
        table_names = []	
        col_names = []
        for table_name, sql in metarec:
            table_names.append(table_name)
            q1 = 'PRAGMA table_info(%s)'
            col_rec = self.get_master_data(con, q1%(table_name))
            for d in col_rec:
                col_names.append(d[1])
            column_names[table_name] = col_names
            col_names=[]
            #print dataset.metadata.column_names[table_name]
	return column_names

    def set_peek( self, dataset, **kwd):
        con = sqlite3.connect(dataset.file_name)
        # no blurb information set for this type
        # rec = self.get_master_data(con, 'select * from metadata where name="blurb"')
        # dataset.blurb = str(dataset.__dict__)
        # dataset.blurb = "<pre>a\r\nb\r\nc</pre>"
        dataset.info = dataset.metadata.config
        rec = self.get_master_data(con, 'select name from sqlite_master where type = "table" and name != "metadata"')
        #rec = self.get_master_data(con, 'select name from sqlite_master where type = "table"')
        table_name = rec[0]
        query = 'select * from %s'
	# limit 7'
        rec = self.get_master_data(con, query%(table_name))
        # get 5 records from the first table in rec
        x=""
        for n in rec:
            x+="%s\n"%(str(n))
        dataset.peek = x
        return True

    # could use display_data like tabular type uses.
    # here is the definition
    # def display_data(self, trans, dataset, preview=False, filename=None, to_ext=None, chunk=None):
    # ...
    #    return trans.fill_template( "/dataset/tabular_chunked.mako",
    #                    dataset = dataset,
    #                    chunk = self.get_chunk(trans, dataset, 0),
    #                    column_number = column_number,
    #                    column_names = column_names,
    #                    column_types = column_types )


class Cm(Db):
    file_ext = "cm"

class Dbi(Db):
    file_ext = "dbi"

    def set_peek( self, dataset, **kwd):
        con = sqlite3.connect(dataset.file_name)
        # no blurb information set for this type
        # rec = self.get_master_data(con, 'select * from metadata where name="blurb"')
        # dataset.blurb = str(dataset.__dict__)
        # dataset.blurb = "<pre>a\r\nb\r\nc</pre>"
        dataset.info = dataset.metadata.config
        #rec = self.get_master_data(con, 'select name from sqlite_master where type = "table" and name != "metadata"')
        rec = self.get_master_data(con, 'select name from sqlite_master where type = "table"')
        table_name = rec[0]
        query = 'select * from %s'
        # limit 7'
        rec = self.get_master_data(con, query%(table_name))
        # get 5 records from the first table in rec
        x=""
        for n in rec:
            x+="%s\n"%(str(n))
        dataset.peek = x
        return True

class Dbm(Db):
    file_ext = "dbm"

