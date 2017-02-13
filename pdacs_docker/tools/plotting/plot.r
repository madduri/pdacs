#!/usr/bin/env Rscript

## This R script plots the following types of graphs
## 1) Plot data out of fof_mf, driver.py script generates, error fractions/error bars
## 2) Plot data out of .cm file, driver.py script generates, error fractions/error plots+line graph
## 3) best fit curve? 
## 4) There might be many more we will extend this tool accordingly or add different R scripts

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RSQLite))


option_list <- list( make_option( c("--verbose", "-v")
                                , action = "store_true"
                                , default = FALSE
                                , help = "Print extra output [%default]"
                                )
                   , make_option( c("--output", "-o")
                                , default = "something.pdf"
                                , help = "Name of PDF output file [%default]"
                                )
                   , make_option( c("--best", "-b")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot best fit curve [%default]"
                                )
                   , make_option( c("--errorbar", "-e")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot error bars [%default]"
                                )
                   , make_option( c("--line", "-l")
                                , action = "store_true"
                                , default = FALSE
                                , help = "plot line curve [%default]"
                                ) 
		   , make_option( c("--table", "-t")
                                , default = "sod_mf"
                                , help = "table from which to get data [%default]"
                                )
                   )

parser <- OptionParser(option_list=option_list, usage="%prog [options] SQLite-file")
args   <- parse_args(parser, positional_arguments = TRUE)
opt    <- args$options

#if (length(args$args) != 2)
#{
#   cat("Incorrect number of required arguments\n\n")
#   print_help(parser)
#   stop()
#}

inputfile  <- as.character(args$args[1])
outputfile <- opt$output

driver <- dbDriver("SQLite")
conn <- dbConnect(driver, dbname=inputfile)
on.exit(dbDisconnect(conn))
on.exit(dbUnloadDriver(drv), add=TRUE)
df <- dbReadTable(conn, opt$table, row.names=NULL)

pdf(file=outputfile)

if(opt$table == 'fof_mf') names(df) <- c('FOF.Mass','clusters','dn.dln','frac.err','iSigM','fsigma','fsigma_fit')
if(opt$table == 'sod_mf') names(df) <- c('SOF.Mass','clusters','cmean','frac.err','variance')
if(opt$table == 'cM') names(df) <- c('SO.Mass', 'clusters','cmean','frac.err','variance')

graph_title <- paste("Plot from ", as.character(opt$table))
 
# Define the top and bottom of the errorbars
tmp.err <- with(df,cmean*frac.err)
limits <- aes(ymax = (cmean + tmp.err), ymin=(cmean - tmp.err))
#legend(1, g_range[2], c("cars","trucks"), cex=0.8, 
 #  col=c("blue","red"), pch=21:22, lty=1:2);

#geom_line(aes(y=variancce/100000), colour='red')
p <- ggplot(df, aes(y=(cmean), x=(SOF.Mass), guide="legend")) + 
  xlab("Mass[Msun/h]") +
  ylab("Concentration") +
  ggtitle(graph_title) + 
  theme(panel.background = element_rect(fill='white', colour='black'), plot.background = element_rect(fill="white")) +
  theme(legend.position = "top")

print(p + geom_point() + scale_x_log10() + geom_point() + geom_errorbar(limits, width=0.025, colour='blue') + geom_line(aes(y=variance/100000), colour='red'))

graphics.off()
