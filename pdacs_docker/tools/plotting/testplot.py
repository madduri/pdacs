#!/usr/bin/env python
#Greg Von Kuster

import sys
from rpy import *

def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()

def column(matrix, i):
    return [row[i] for row in matrix]

def main():

    in_fname = sys.argv[1]
    out_fname = sys.argv[2]
    print in_fname
    print out_fname
    title = sys.argv[3]
    xlab = sys.argv[4]
    ylab = sys.argv[5]

    matrix = [[0.0255, 0.4462], [0.0302, 0.4713], [0.034, 0.4984], [0.041, 0.5278], [0.0477, 0.5599], [0.0568, 0.5958], [0.0621, 0.636], [0.0745, 0.6771], [0.0867, 0.7236], [0.1015, 0.7762], [0.127, 0.8321], [0.169, 0.8996], [0.1826, 0.9592], [0.2673, 1.0529], [0.3333, 1.1533], [0.3333, 1.2389]]
    print column(matrix, 0)
    print column(matrix, 1)
   
    r.pdf( out_fname, 8, 8 )
    print array(matrix)
    r.plot(column(matrix, 0), column(matrix, 1), type="p", main=title, xlab=xlab, ylab=ylab, col="blue", pch=19 )
    print "didnt plot"
    r.dev_off()
    r.quit( save="no" )

if __name__ == "__main__":
    main()
