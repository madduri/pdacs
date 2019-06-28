#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RSQLite))

args = commandArgs()
print(args)


#driver = dbDriver("SQLite")
#conn = dbConnect(driver, dbname="${out_file1}")
#on.exit(dbDisconnect(conn))
#on.exit(dbUnloadDriver(drv), add=TRUE)
#ires = dbSendQuery(conn=conn, "create table metadata(key TEXT, value NUMERIC)")
#dbClearResult(ires)
      #for $i, $s in enumerate( $series )
#          s${i} = read.table( "${s.input.file_name}", sep="\t", header=TRUE)
#          tempname = "${s.tabname}"
#          print(tempname)
#          if (nchar(tempname) == 0) {
#            tabname1 = strsplit("${s.input.name}", " on")
#            tabname2 = strsplit("${s.input.file_name}", "_")
#            tempname = paste(tabname1[[1]][1], tabname2[[1]][2])
#          }
#          tablename = toString(strsplit(gsub('([[:punct:]])|\\s+','_', tempname), ".dat"))
#          ires = dbWriteTable(conn=conn, name=tablename, value=as.data.frame(s${i}), row.names=FALSE, header=TRUE)
      #end for





