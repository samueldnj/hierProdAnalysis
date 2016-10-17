# --------------------------------------------------------------------------
# tools.R
# 
# Script that contains tool functions for the coastwide multispecies 
# simulation-estimation procedure.
# 
# Author: Samuel D N Johnson
# Date: 24 August, 2016
#
#
#	List of tools:
#		lisread()			Function that reads ADMB style files as a list object
#   writeADMB()   Function that writes ADMB dat/pin files from list objects
#   saveSim()     Function that saves the output from a simulation run
# 
# --------------------------------------------------------------------------

saveSim <- function(blob, folder, path)
{
  rDataFile <- paste ( folder, ".RData", sep = "" )
  dataFilePath <- file.path ( path, rDataFile )
  save ( blob, file = dataFilePath )
  cat( "\nMSG (saveSim) Saved ",rDataFile,"in project/",folder,"/.\n" )
}


# lisread()
# lisread: Function to read a list of data objects from a file.
# The initial characters "##" denote a comment line (ignored).
# The initial characters "# " denote a variable name.
# All other lines must contain scalars or vectors of numbers.
# Furthermore, all rows in a matrix must contain the same number of
# columns. Row and column vectors are not converted to matrices.
#
# fname  : File name.
# quiet  : If true, shut up about reporting progress.
# result : List object with components named in the file.

# Original functions courtesy of Jon Schnute.
# Modifications by A.R. Kronlund.
lisread <- function( fname,quiet=TRUE )
{
  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x,white=T )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )
    xc
  }

  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables

  if( nvar==1 )
    nrow <- c( nf2 + 1 ) - ilab - 1
  else
    nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1     # number of rows in var i
  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1
    i2 <- i1 + nrow[i] - 1                       # range of lines for var i
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}


# Write element i of list object x to the active file. Note that any
# elements must be numeric (as this is the only thing ADMB can read in),
# also sub-list structure is removed and elements of lists are written as 
# vectors, probably to be read in to ADMB as ragged arrays
# usage:
# cat ( "## writeAMDB output file, created ", format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", sep = "", file = activeFile )
# lapply ( X = seq_along(x), FUN = writeADMB, x, activeFile)
writeADMB <- function( i, x, activeFile ) {
  # Write element name
  cat( paste( "# ", names(x)[[i]], "\n", sep = "" ), file=activeFile, append=TRUE )
  # Write element value
  if( is.matrix(x[[i]]) ) {
    # Matrices
    write.table( x[[i]], file=activeFile, row.names=FALSE, col.names=FALSE, append=TRUE )
  } else if( is.list(x[[i]]) ) {
    # Lists
    lapply( x[[i]], "\n", FUN=cat, file=activeFile, append=TRUE )
  } else {
    # Vectors/Scalars
    cat( x[[i]], "\n", file=activeFile, append=TRUE )  
  }
}  # end writeADMB

# Read ragged arrays from ADMB rep object
readRagged <- function( x, var, ncols ) {
  outList <- list()
  for( i in 1:length(ncols) ) {
    outList[[i]] <- numeric(ncols[i])
    for( j in 1:ncols[i] ) {
      varName <- paste( var, "_", i, j, sep=""  )
      k <- which( names(x)==varName )
      outList[[i]][j] <- x[[k]]
    }
  }
  outList
}   # end readRagged()

read.fit <-
function(ifile)
{
  # __Example:             
  # file <-("~/admb/simple")
  # A <- reptoRlist(file)
  # Note there is no extension on the file name.
  
  ## The following is a contribution from:
  ## Anders Nielsen that reads the par & cor files.
  ret<-list() 
  parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
   what='', n=16, quiet=TRUE)[c(6,11,16)]) 
  ret$nopar<-as.integer(parfile[1]) 
  ret$nlogl<-parfile[2] 
  ret$maxgrad<-parfile[3] 
  file<-paste(ifile,'.cor', sep='') 
  lin<-readLines(file) 
  ret$npar<-length(lin)-2 
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
  ret$names<-unlist(lapply(sublin,function(x)x[2])) 
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
  ret$cor<-matrix(NA, ret$npar, ret$npar) 
  corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
  ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
  ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  return(ret)
}