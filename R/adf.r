#' Create an abstract data frame
#'
#' @importFrom           iotools mstrsplit dstrsplit chunk.reader chunk.apply
#' @param description    character string.  A description of the connection; the
#                          path to a file for file connections.
#' @param conMethod      string indicating the connection method.
#' @param sep            character seperating the data columns. Ignored when
#'                         chunkFormatter is given.
#' @param nsep           character seperating the column of rownames. Set to
#'                         NA to generate automatic rownames. Ignored when
#'                         chunkFormatter is given.
#' @param strict         logical. Whether the parser should run in strict mode.
#'                         Ignored when chunkFormatter is given.
#' @param header         logical indicating whether the first line of data, after
#'                         skip if provided, contains variable names. Ignored if
#'                         chunkFormatter is also provided, as adf would not know how
#'                         to parse this first row.
#' @param colNames       an optional character vector of column names. If missing
#'                         and header=TRUE, these will be determined be the first
#'                         row of data. Otherwise, names will be constructed by
#'                         pasting the character 'V' with the column number.
#' @param colClasses     an optional character vector of column classes. If
#'                         named and colNames is missing, the names will be
#'                         used for colNames. If missing, will be automatically
#'                         determined.
#' @param levels         a named list, with names coorisponding to the colNames
#'                         of class character (or factor). Each element gives the
#'                         levels for the coorisponding variable. These will be
#'                         automatically determined when missing.
#' @param skip           number of lines to strip off of the connection before
#'                         parsing
#' @param nrowsClasses   number of rows to pull for determining colClasses and
#'                         factor levels when needed
#' @param chunkProcessor a function to apply post-processing to a formatted
#'                         dataframe. Usually only used without a chunkFormatter.
#' @param chunkFormatter an optional function turning the raw connection into
#'                         a dataframe. It must accept four parameters: data,
#'                         colNames, colClasses, levels. If missing this will
#'                         be constructed automatically.
#' @param minSplits an optional paramater for use with SparkR.
#'
#' @return An abstract data frame object, which can be used with iolm, ioglm, and
#'         other related functions.
#' @export
adf = function(description, conMethod = c("file", "gzfile", "bzfile", "xzfile"),
                sep="|", nsep=NA, strict=TRUE, header=FALSE,
                colNames, colClasses, levels = list(),
                skip=0L, nrowsClasses=100L, chunkProcessor = identity,
                chunkFormatter, minSplits=2L) {

  if (missing(colNames) && !missing(colClasses) && !is.null(names(colClasses)))
    colNames = names(colClasses)
  if (header && !missing(chunkFormatter)) {
    warning("Header not used with user-defined chunkFormatter")
    header = FALSE
  }
  output = list(chunkProcessor = chunkProcessor, skip = skip + header, levels=list())

  if (!missing(colNames))
    output$colNames = colNames
  if (!missing(colClasses))
    output$colClasses = colClasses

  # Construct the connection method
  if (is.function(conMethod)) {
    output$createNewConnection = conMethod
  } else if (inherits(conMethod,"jobj")) {
    # output$createNewConnection = SparkR::textFile(conMethod, description, minSplits)
  } else {
    description = path.expand(description)
    conMethod = match.arg(conMethod)
    output$createNewConnection = call(conMethod, description=description, open="rb")
  }

  # Build the chunkFormatter
  needHead = (missing(colNames) || missing(colClasses) ||
              length(setdiff(colNames[colClasses %in% c("character", "factor")],names(levels))))
  if (needHead) {
    if (inherits(output$createNewConnection, "RDD")) {
      # rows = SparkR::lapplyPartition(output$createNewConnection, function(z) head(z,n=nrowsClasses))
      # rows = SparkR::collect(rows)
      # rows = as.vector(rows, mode="character")
    } else {
      on.exit(close(con))
      if (is.function(output$createNewConnection)) {
        con = output$createNewConnection()
      } else {
        con = eval(output$createNewConnection)
      }
      readLines(con, n=output$skip - header)
      if (missing(colNames) && header)
        output$colNames = strsplit(readLines(con,1L), split=sep, fixed=TRUE)[[1]]

      rows = readLines(con, n=nrowsClasses)
    }

    if (missing(chunkFormatter)) {
      z = as.data.frame(iotools::mstrsplit(rows, sep=sep, nsep=nsep),
                        stringsAsFactors=FALSE)
      if (missing(colNames) && !header)
        output$colNames = sprintf("V%d", seq(ncol(z)))
    } else {
      z = chunkFormatter(rows,
                          ifelse(missing(colNames), NULL, colNames),
                          ifelse(missing(colClasses), NULL, colClasses),
                          output$levels)
      output$colNames = names(z)
    }
    if (missing(colClasses))
      output$colClasses = sapply(z,function(v) class(type.convert(v, as.is=TRUE)))
    output$colClasses[output$colClasses == "factor"] = "character"
    for (j in which(output$colClasses == "character")) {
      if (is.null(output$levels[[output$colNames[j]]]))
        output$levels[[output$colNames[j]]] = levels(factor(z[,j]))
    }
  }

  if (length(output$colNames) != length(output$colClasses))
    stop(sprintf("Number of column names (%d) not equal to number of column classes (%d)",
         length(output$colNames), length(output$colClasses)))

  if (!missing(chunkFormatter)) {
    output$chunkFormatter = chunkFormatter
  } else {
    output$chunkFormatter = default.chunkFormatter(sep=sep, nsep=nsep, strict=strict)
  }

  class(output) = c("adf")
  return(output)
}

#' Construct a default chunkFormatter
#'
#' @param sep     character seperating the data columns.
#' @param nsep    character seperating the column of rownames. Set to
#'                  NA to generate automatic rownames.
#' @param strict  logical. Whether the parser should run in strict mode.
default.chunkFormatter = function(sep=sep, nsep=nsep, strict=strict) {
  function(data, colNames, colClasses, levels) {
    col_types = colClasses
    names(col_types) = colNames
    if (is.list(data))
      data = as.vector(data, mode="character")
    z = dstrsplit(data, col_types=col_types, sep=sep, nsep=nsep, strict=strict)
    for (j in which(colClasses == "character")) {
      if (!is.null(levels[[colNames[j]]]))
        z[,j] = factor(z[,j], levels=levels[[colNames[j]]])
    }
    return(z)
  }
}

#' Update abstract data frame with all factor levels
#'
#' @param x   an abstract dataframe
#'
#' @return a new abstract data frame, with (possibly) new
#'         factor levels derived from scanning the entire
#'         dataset.
#' @export
allFactorLevels = function(x) {
  if (!inherits(x, "adf")) stop("x must be an 'adf' object!")

  index = which(x$colClasses == "character")
  if (length(index) == 0L) return(x)

  x$levels = vector(mode="list", length=length(x$colClasses))
  names(x$levels) = x$colNames
  z = adf.apply(x, FUN = function(d,passedVars) {
      lapply(d[,index],unique)
    })
  for (j in index) {
    x$levels[[x$colNames[j]]] =
      sort(unique(unlist(Map(function(w) w[x$colNames[j]], z),
                                    use.names=FALSE)))
  }
  return(x)
}

#' Convert a data.frame to an abstract data frame
#'
#' @param x    a data.frame
#' @param ...  other arguments. Currently unused
#'
#' @return an abstract data frame object
#' @export
as.adf = function(x, ...) {
  output = list()

  output$levels = lapply(x, levels)
  factor_cols = which(sapply(x, class) == "factor")
  for (fc in factor_cols) x[,fc] = as.character(x[,fc])

  output$createNewConnection = iotools::as.output.data.frame(x)
  output$skip = 0L
  output$chunkProcessor = identity
  output$chunkFormatter = default.chunkFormatter("|", NA, TRUE)
  output$colNames = names(x)
  output$colClasses = sapply(x, class)
  output$levels = lapply(x, levels)
  for (j in setdiff(which(output$colClasses == "character"), factor_cols))
    output$levels[[j]] = unique(x[,j])

  class(output) = c("adf")
  return(output)
}

#' Low level function for applying over an abstract data frame
#'
#' @param x                an abstract data.frame
#' @param FUN              function to apply over each chunk; must accept two
#'                         inputs: the data.frame (or model list), and the
#'                         input to passedVars
#' @param type             type of data to give as an input to FUN. If model
#'                         or sparse model, this is a list giving the response (y),
#'                         model matrix (x), weights (w), and offset (offset)
#'                         from the input forumal.
#' @param formula          a formula to use with type equal to model or sparse.model
#' @param contrasts        contrasts to use with type equal to model or sparse.model
#' @param subset           a string to to use with type equal to model or sparse.model.
#'                         Will be evaluated in the environment of the data frame
#'                         (ex. subset="V2 + V3 > V4")
#' @param weights          a string to to use with type equal to model or sparse.model.
#'                         Will be evaluated in the environment of the data frame.
#' @param na.action        a function which indicates what should happen when the data
#'                         contain 'NA's. See lm.fit for more details.
#' @param offset           a string to to use with type equal to model or sparse.model.
#'                         Will be evaluated in the environment of the data frame.
#' @param passedVars       Option list of arguments which is passed to FUN
#' @param ...              Other arguments to pass.
#' @param chunk.max.line   integer. Maximum number of lines to parse in a single chunk.
#' @param CH.MAX.SIZE      integer. Maximum size in bytes of a chunk.
#' @param CH.MERGE         method for combining the output of chunks. Typically list or rbind.
#' @param parallel         integer. the number of parallel processes to use in the calculation
#'                         (*nix only).
#'
#' @return an object of class type
#' @export
adf.apply = function(x, FUN, type=c("data.frame", "model", "sparse.model"),
                     formula=NULL, contrasts=NULL, subset=NULL, weights=NULL,
                     na.action=NULL, offset=NULL, passedVars=NULL, ...,
                     chunk.max.line=65536L, CH.MAX.SIZE=33554432L,
                     CH.MERGE = list, parallel=1L) {

  if (!inherits(x, "adf")) stop("x must be an 'adf' object!")
  type = match.arg(type)

  FUN2 = function(z) {
    df = x$chunkFormatter(z, x$colNames, x$colClasses, x$levels)
    df = x$chunkProcessor(df)
    if (type == "data.frame") return(FUN(df,passedVars))

    if (!is.null(subset)) {
      subset = eval(parse(text=paste0("with(df, ", subset ,")")))
      df = df[subset,]
    }

    IO_OFFSET = IO_WEIGHTS = NULL # to silence NOTES from R CMD CHECK
    environment(formula) = parent.frame()
    if (is.null(weights) && is.null(offset)) {
      mf = model.frame(formula=formula, data=df)
    } else if (is.null(weights)) {
      df$IO_OFFSET = eval(parse(text=paste0("with(df, ", offset ,")")))
      mf = model.frame(formula=formula, data=df, offset=IO_OFFSET)
    } else if (is.null(offset)) {
      df$IO_WEIGHTS = eval(parse(text=paste0("with(df, ", weights ,")")))
      mf = model.frame(formula=formula, data=df, weights=IO_WEIGHTS)
    } else {
      df$IO_OFFSET = eval(parse(text=paste0("with(df, ", offset ,")")))
      df$IO_WEIGHTS = eval(parse(text=paste0("with(df, ", weights ,")")))
      mf = model.frame(formula=formula, data=df, offset=IO_OFFSET, weights=IO_WEIGHTS)
    }
    mt = attr(mf, "terms")

    y = model.response(mf, "numeric")
    w = as.vector(model.weights(mf))
    offset = as.vector(model.offset(mf))
    if (type == "sparse.model")
      x = Matrix::sparse.model.matrix(mt, mf, contrasts.arg=contrasts)
    else
      x = model.matrix(mt, mf, contrasts.arg=contrasts)

    FUN(list(y=y,x=x,w=w,offset=offset), passedVars)
  }

  if (inherits(x$createNewConnection, "RDD")) {
    # output = SparkR::lapplyPartition(x$createNewConnection, function(z) list(FUN2(z)))
    # output = SparkR::collect(output)
  } else {
    on.exit(close(con))
    if (is.function(x$createNewConnection)) {
      con = x$createNewConnection()
    } else if (is.raw(x$createNewConnection)) {
      con = rawConnection(x$createNewConnection)
    } else {
      con = eval(x$createNewConnection)
    }
    if (x$skip > 0L) readLines(con, n=x$skip)
    cr = chunk.reader(con, max.line=chunk.max.line)

    output = chunk.apply(cr, FUN2, CH.MERGE = CH.MERGE, CH.MAX.SIZE = CH.MAX.SIZE,
                         parallel=parallel)
  }

  return(output)
}

#' Construct model.frame
#'
#' A custom version of lm to construct a model frame that
#' does not drop unused levels (because these should be)
#' in other chunks.
#'
#' @param formula          a formula to use with type equal to model or sparse.model
#' @param data             a data.frame object
#' @param subset           a string to to use with type equal to model or sparse.model.
#'                         Will be evaluated in the environment of the data frame
#'                         (ex. subset="V2 + V3 > V4")
#' @param weights          a string to to use with type equal to model or sparse.model.
#'                         Will be evaluated in the environment of the data frame.
#' @param na.action        a function which indicates what should happen when the data
#'                         contain 'NA's. See lm.fit for more details.
#' @param model            an optional subset vector
#' @param contrasts        contrasts to use with type equal to model or sparse.model
#' @param offset           an optional offset
lm.model.frame =
function (formula, data, subset, weights, na.action,
    contrasts = NULL, offset)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
      "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- FALSE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  return(mf)
}

#' Print the details of an abstract data frame
#'
#' @method print adf
#' @param x    a data.frame
#' @param ...  other arguments. Currently unused
#'
#' @export
print.adf = function(x, ...) {
  w = min(max(nchar(x$colNames), 10L), 40L)
  cat(sprintf("    An abstract data frame with %d columns:\n\n", length(x$colClasses)))
  cat(sprintf(paste0("    %-", w, "s  %-10s"), x$colNames, x$colClasses),sep="\n")
}
