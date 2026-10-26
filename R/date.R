## not used.
## we don't guess mrsd from tip labels
## user should provide mrsd parameter in ggtree if they want to use time-scaled tree
scaleX_by_time <- function(df, as.Date=FALSE) {
    time <- as.numeric(with(df, gsub(".*[_/]{1}(\\d+\\.*\\d+)$", "\\1", label[isTip])))
    latest <- which.max(time)

    scaleX_by_time_from_mrsd(df, decimal2Date(time[latest]), as.Date)
}

scaleX_by_time_from_mrsd <- function(df, mrsd, as.Date) {
    if (!inherits(mrsd, "Date")) mrsd <- as.Date(mrsd)
    date <- Date2decimal(mrsd)

    df$x <- df$x + date - max(df$x)
    df$branch <- with(df, (x[match(parent, node)] + x)/2)
    ## df$branch <- (df[df$parent, "x"] + df[, "x"])/2

    if (as.Date) {
        df$x <- decimal2Date(df$x)
        df$branch <- decimal2Date(df$branch)
    }

    return(df)
}


##' convert Date to decimal format, eg "2014-05-05" to "2014.34"
##'
##' @param x a Date, or a numeric that is passed through unchanged.
##' @return numeric year with fraction, e.g. 2014.34.
##' @examples
##' Date2decimal(as.Date("2014-05-05"))
##' @export
##' @author Guangchuang Yu
Date2decimal <- function(x) {
    if (is(x, "numeric")) {
        return(x)
    }

    if (is(x, "character")) {
        x <- as.Date(x)
    }
    year <- format(x, "%Y")
    y <- x - as.Date(paste0(year, "-01-01"))
    year <- as.numeric(year)
    ndate <- ifelse(year %% 4 == 0, 366, 365)
    year + as.numeric(y)/ndate
}

##' convert decimal format to Date, eg "2014.34" to "2014-05-05"
##'
##' @param x numeric year with fraction, e.g. 2014.34.
##' @return a Date.
##' @examples
##' decimal2Date(2014.34)
##' @export
##' @author Guangchuang Yu
decimal2Date <- function(x) {
    date <- as.Date(paste0(floor(x), "-01-01"))
    date + as.numeric(sub("^\\d+", "0", x)) * 365
}


