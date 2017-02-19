# binsize
min_binsize = 10

# lower boundary
lowerboundary <- function(x, increment) {
    return (floor(x / increment) * increment)
}

# upper boundary
upperboundary <- function(x, increment) {
    return (ceiling(x / increment) * increment)
}

# round to decimals
roundup <- function(x) {
    return (sign(x) * 10^ceiling(log10(abs(x))))
}

# wrapper
wrapper <- function(table, columns, options) {

    # get binsize
    binsize = max(as.integer(options$binsize), min_binsize)
    
    # initialize output list
    l <- list()

    # loop through all columns
    m <- list()
    for (key in names(columns)) {
        # load column data
        column <- as.numeric(columns[key])
        column_data <- sapply( table[column], as.numeric )

        # collect vectors in list
        m <- append(m, list(column_data))
    }
    
    # get min/max boundaries
    min_value <- min(unlist(m))
    max_value <- max(unlist(m))
    
    # identify range
    diff <- max_value - min_value
    
    # identify increment
    increment <- roundup(diff / binsize)
    
    # fix min value
    min_value <- lowerboundary(min_value, increment)
    max_value <- upperboundary(max_value, increment)
    
    # update range
    diff <- max_value - min_value
    
    # fix bin size
    binsize = round(diff / increment)
    
    # fix max value
    max_value <- min_value + binsize * increment
    
    # check if single bin is enough
    if (min_value == max_value) {
        l <- append(l, max_value)
        for (key in seq(m)) {
            l <- append(l, 1.0)
        }
        return (l)
    }
    
    # fix range and bins
    bin_seq = seq(min_value, max_value, by=increment)
    
    # add as first column
    l <- append(l, list(bin_seq[2: length(bin_seq)]))
    
    # loop through all columns
    for (key in seq(m)) {
        # load column data
        column_data <- m[[key]]
        
        # create hist data
        hist_data <- hist(column_data, breaks=bin_seq, plot=FALSE)
        
        # normalize densities
        count_sum <- sum(hist_data$counts)
        if (count_sum > 0) {
            hist_data$counts = hist_data$counts / count_sum
        }

        # collect vectors in list
        l <- append(l, list(hist_data$counts))
    }
    
    # return
    return (l)
}
