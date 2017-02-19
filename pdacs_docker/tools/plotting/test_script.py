#!/usr/bin/env python

import sqlite3 as lite
con1 = lite.connect("/global/project/projectdirs/hacc/PDACS/working/interim.db")
cur1 = con1.cursor()
cur1.execute("drop table if exists forplot")
cur1.execute("create table forplot(x_col NUMERIC, y_col NUMERIC)")
#for $i, $s in enumerate( $series )
##if ${s.format_type}=="cM":
con = lite.connect("${s.input.file_name}")
##if ${s.format_type}=="kpk":
##  con = lite.connect("${s.kpk_input.file_name}")
##if ${s.format_type}=="db":
##  con = lite.connect("${s.db_input.file_name}")
cur = con.cursor()
cur.execute("select ${s.xcol} from ${s.table_name}")
x=cur.fetchall()
cur.execute("select ${s.ycol} from ${s.table_name}")
y = cur.fetchall()
for d,e in zip(x, y):
   int1 = d[0]
   int2 = e[0]
   print int1, int2
print "Table 1 done and added to combined file"
##for a, b in zip(x${i}, y${i}):
##  tuples = (a, b)
##  cur1.execute("insert into forplot values(?,?)", tuples)
#end for
##con1.commit()

