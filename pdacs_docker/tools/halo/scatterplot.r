#!/usr/bin/env Rscript

## Read the output of cM, in its SQLite database format, and make a plot of log(Mass) vs. concentration.
##

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(RSQLite))

option_list <- list( make_option( c("--verbose", "-v")
                                , action = "store_true"
                                , default = FALSE
                                , help = "Print extra output [%default]"
                                )
                   , make_option( c("--output", "-o")
                                , default = "scatterplot.pdf"
                                , help = "Name of PDF output file [%default]"
                                )
                   , make_option( c("--formula", "-f")
                                , default = "y ~ x"
                                , help = "Formula used in the xyplot command [%default]"
                                )
                   , make_option( c("--loess", "-l")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot loess curve [%default]"
                                )
                   , make_option( c("--regression", "-r")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot linear regression [%default]"
                                )
                   , make_option( c("--nopoints", "-n")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot data points curve [%default]"
                                )
                   , make_option( c("--table", "-t")
                                , default = "data"
                                , help = "table from which to get data [%default]"
                                )
                   )

parser <- OptionParser(option_list=option_list, usage="%prog [options] SQLite-file")
args   <- parse_args(parser, positional_arguments = TRUE)
opt    <- args$options

if (length(args$args) != 1)
{
   cat("Incorrect number of required arguments\n\n")
   print_help(parser)
   stop()
}

inputfile  <- as.character(args$args)
outputfile <- opt$output

if (opt$verbose)
{
  cat("Will process ", inputfile, " and write to ", outputfile, "\n")
  cat("Will plot formula ", opt$formula, "\n")  
  cat("Will get data from ", opt$table, "\n")
  if (opt$loess)      cat("Will plot loess curve\n")
  if (opt$regression) cat("Will plot linear regression\n")
  if (opt$nopoints)   cat("Will plot data points\n")
}

driver <- dbDriver("SQLite")
conn <- dbConnect(driver, dbname=inputfile)
on.exit(dbDisconnect(conn))
on.exit(dbUnloadDriver(drv), add=TRUE)

df <- dbReadTable(conn, opt$table, row.names=NULL)
if (opt$nopoints)   type.plot <- NULL else type.plot <- c("p")
if (opt$loess)      type.plot <- append(type.plot, "smooth")
if (opt$regression) type.plot <- append(type.plot, "r")

pdf(file=outputfile)
xyplot( as.formula(opt$formula),
        df,
        type = type.plot,
        col.line = "red",
        grid=TRUE)
graphics.off()







