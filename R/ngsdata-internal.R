
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Tue Feb 09 13:09:40 CST 2021 -0600 (Week 06)
##
## Contents:
##   Helper functions to create a sample-by-feature matrix
## 
## Reference: 
## 
## ************************************************************************

stop.not.in.vector <- function(
  x,
  y
  ){
  .x <- deparse(substitute(x))
  .y <- deparse(substitute(y))
  .xv <- vconcat(x)
  .yv <- vconcat(y)
  if (!(x %in% y)) {
    .msg <- lprintf("`%(.x)s` [%(.xv)s] must be in vector `%(.y)s` [%(.yv)s].")
    stop(.msg)
  }
}



get_toy__data.score <- function(
  ){
  system.file("extdata","data.sample.txt",package="NGSdata",mustWork=TRUE)
}

get_toy__data.region <- function(
  ){
  system.file("extdata","data.feature.txt",package="NGSdata",mustWork=TRUE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @return 
##' @author
##' @examples
##' get_FUNC.SCORE_vconcat(c(1:6,99))
##' 
get_FUNC.SCORE_vconcat <- function(x){
  vconcat(x,capsule=FALSE,quote=FALSE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @return 
##' @author 
##' get_FUNC.SCORE_mean(c(1:6,99))
get_FUNC.SCORE_mean <- function(x){
  .func <- function(x){mean(x,na.rm=TRUE)}
  .func(x)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @return 
##' @author 
##' get_FUNC.SCORE_median(c(1:6,99))
get_FUNC.SCORE_median <- function(x){
  .func <- function(x){median(x,na.rm=TRUE)}
  .func(x)
}

  

##' Read a file in `data.table` format
##'
##' 
##' @title Read data as a `data.table` 
##' @param inFpath file, See \code{read.table}
##' @param header See \code{read.table}
##' @param ... 
##' @return 
##' @author
##' @examples
##' 
read.data <- function(
  inFpath,
  header=TRUE,
  ...
  ){
  ret <- read.table(inFpath,stringsAsFactors=FALSE,header=header,...)
  setDT(ret)
  return(ret)
}

data.as.matrix <- function(
  x,
  formula=sampleName ~ regionName,
  which.col.score = "score"
  ){
  ret <- reshape2::acast(x,formula,value.var=which.col.score)
  return(ret)
}


##' `data.region` is a data table of regions to score.
##'
##' e.g. a data table of region information by feature (e.g. gene).
##' @title 
##' @param x 
##' @param which.col.chr 
##' @param which.col.start 
##' @param which.col.end 
##' @param which.col.name 
##' @return 
##' @author 
prep_data.region <- function(
  x,
  which.col.chr = "chr", 
  which.col.start = "start", 
  which.col.end = "end",
  which.col.regionName = "name"
  ){
  .colnames <- c(which.col.chr,which.col.start,which.col.end,which.col.regionName)
  if (!all(.colnames %in% colnames(x))){
    sapply(.colnames,stop.not.in.vector,y=colnames(x))
  }
  x <- x[,.colnames,with=FALSE]
  setnames(x,.colnames,c("chr","start","end","regionName"))
  x[,chr:=as.character(chr)]
  x[,start:=as.numeric(start)]
  x[,end:=as.numeric(end)]
  x[,regionName:=as.character(regionName)]
  return(x)
}

##' `data.score` is a data table of score information over regions
##'
##' e.g. a data table of score information over regions by sample.
##' @title To prepare `data.score`
##' @param x 
##' @param which.col.chr colname for chromosomes or other sequences
##' @param which.col.start colname for the start of a region
##' @param which.col.end colname for the end of a region
##' @param which.col.name colname for the name of a region (e.g. sample names)
##' @param which.col.score colname for the score of a region
##' @return 
##' @author 
prep_data.score <- function(
  x,
  which.col.chr = "chr", 
  which.col.start = "start", 
  which.col.end = "end",
  which.col.sampleName = "name",
  which.col.score = "score"
  ){
  .colnames <- c(which.col.chr,which.col.start,which.col.end,which.col.sampleName,which.col.score)
  if (!all(.colnames %in% colnames(x))){
    sapply(.colnames,stop.not.in.vector,y=colnames(x))
  }
  x <- x[,.colnames,with=FALSE]
  setnames(x,.colnames,c("chr","start","end","sampleName","score"))
  x[,chr:=as.character(chr)]
  x[,start:=as.numeric(start)]
  x[,end:=as.numeric(end)]
  x[,sampleName:=as.character(sampleName)]
  x[,score:=as.numeric(score)]
  return(x)
}
