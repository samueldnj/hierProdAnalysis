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
# ToDo: 1. Complete list of tools
#       2. Write headers for saveSim()
#
#	List of tools:
#		lisread()			Function that reads ADMB style files as a list object (JS)
#   writeADMB()   Function that writes ADMB dat/pin files from list objects
#   saveSim()     Function that saves the output from a simulation run
#   read.admb()   Function that reads ADMB output (SJDM)
# 
# --------------------------------------------------------------------------

# panLegend   (Place legend in plot region)
# Purpose:    Place a legend in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "legend" to be
#             passed.
# Parameters: x, y are the coordinates of the legend
#             legTxt is the text associated with the legend
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLegend <- function( x, y, legTxt, ... )
{
  # Allows legend to be placed at 0<x<1, 0<y<1.
  usr  <- par( "usr" )
  yLog <- par("ylog")
  xLog <- par("xlog")
  # Check for log-transformed axes and adjust usr commands as needed
    # note: when a log scale is in use, 
    #           usr gives limits in the form 10 ^ par("usr")
  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE) {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
     usr[3:4]<-10 ^ par("usr")[3:4]
     par( usr=c(0,1,0,1), ylog=FALSE )
   } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
     usr[1:2]<-10 ^ par("usr")[1:2]
     par( usr=c(0,1,0,1), xlog=FALSE )
   } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
     usr[1:4]<-10 ^ par("usr")[1:4]
     par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
   } 
  legend( x, y, legend=legTxt, ... )
  par( usr=usr )
  return( NULL )
}

# panLab      (Place text labels in plot region)
# Purpose:    Place a text label in the plot region defined by (0,1), (0,1).
#             The ... notation allows all parameters available to "text" to be
#             passed.
# Parameters: x, y are the coordinates of the label
#             txt is the text
# Returns:    NULL (invisibly)
# Source:     A.R. Kronlund
# Revised:    K.Holt; 13-Jan-10 to accomodate axes on log scale
panLab <- function( x, y, txt, ... )
{
  # Allows text to be placed in plot panel at 0<x<1, 0<y<1.
  usr <- par( "usr" )
  
  yLog <- par("ylog")
  xLog <- par("xlog")
  
  # Check for log-transformed axes and adjust usr commands as needed
  # note: when a log scale is in use, 
  #           usr gives limits in the form 10 ^ par("usr")

  # Case 1: neither axis is on the log scale
  if (yLog==FALSE & xLog==FALSE)
  {
    par( usr=c(0,1,0,1) )
  }
  # Case 2: only the y-axis is on log scale
  if (yLog==TRUE & xLog==FALSE) 
  {
    usr[3:4]<-10 ^ par("usr")[3:4]
    par( usr=c(0,1,0,1), ylog=FALSE )
  } 
  # Case 3: only the x-axis is on log scale
  if (yLog==FALSE & yLog==TRUE) 
  {
    usr[1:2]<-10 ^ par("usr")[1:2]
    par( usr=c(0,1,0,1), xlog=FALSE )
  } 
  # Case 4: both axes are on the log scale
  if (yLog==TRUE & xLog==TRUE) 
  {
    usr[1:4]<-10 ^ par("usr")[1:4]
    par( usr=c(0,1,0,1), xlog=FALSE, ylog=FALSE )
  } 
  text( x, y, txt, ... )
  par( usr=usr )
  return( NULL )
}

# loadSim()
# Loads the nominated simulation blob into memory, so that plot functions
# run faster.
# inputs:   sim=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
loadSim <- function (sim=1)
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path="./project",full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  simList <- dirList[grep(pattern="sim",x=dirList)]
  folder <- simList[sim]

  # Load the nominated blob
  blobFileName <- paste(folder,".RData",sep="")
  blobPath <- file.path(getwd(),"project",folder,blobFileName)
  load ( file = blobPath )

  assign( "blob",blob,pos=1 )
  cat("MSG (loadSim) Simulation ", folder, " loaded from ./project/\n", sep="" )
}


# calcDetFit()
# For troubleshooting bad model fits, by stepping through the
# empirical covariance matrix and calculating determinants for sub-matrices.
# Hopefully, the sub-matrices without full rank can be identified and 
# the culprits can be found.
# inputs:   fit=list returned by read.fit
# outputs:  det=named vector of deteterminants for submatrices missing
#               the named variable
calcDetFit <- function (fit)
{
  npar  <- fit$npar
  cor   <- fit$cor

  subDet <- numeric ( length = npar )

  for(j in 1:npar)
  {
    subMat <- cor[-j,-j]
    subDet[j] <- det(subMat)
  }
  return(subDet)
}

# saveSim()
# Function that....
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
read.admb <-
function(ifile)
{ 
  ret=read.fit(ifile)
  
  fn=paste(ifile,'.rep', sep='')
  A=read.rep(fn)
  A$fit=ret
  
  pfn=paste(ifile,'.psv',sep='')
  if(file.exists(pfn))
    A$post.samp=read.psv(pfn)
  
  return(A)
}

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

read.rep <- 
function(fn)
{
  # The following reads a report file
  # Then the 'A' object contains a list structure
  # with all the elemements in the report file.
  # In the REPORT_SECTION of the AMDB template use 
  # the following format to output objects:
  #   report<<"object \n"<<object<<endl;
  #
  # The part in quotations becomes the list name.
  # Created By Steven Martell
  options(warn=-1)  #Suppress the NA message in the coercion to double
  
  
  ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))

    if(is.numeric(dum))#Logical test to ensure dealing with numbers
    {
      A[[vnam[i]]]=dum
    }
  }
  options(warn=0)
  
  return(A)
}

read.psv <-
function(fn, nsamples=10000)
{
  #This function reads the binary output from ADMB
  #-mcsave command line option.
  #fn = paste(ifile,'.psv',sep='')
  filen <- file(fn, "rb")
  nopar <- readBin(filen, what = integer())
  mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
  mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
  close(filen)
  return(mcmc)
}

gletter <-
function(i=1)
{
  usr=par("usr"); inset.x=0.05*(usr[2]-usr[1]); inset.y=0.05*(usr[4]-usr[3])
  text(usr[1]+inset.x,usr[4]-inset.y,paste("(",letters[i],")",sep=""),cex=1.,font=1)
}