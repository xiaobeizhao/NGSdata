
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Tue Feb 09 13:36:53 CST 2021 -0600 (Week 06)
## 
## 
## Reference: 
## Xmisc/Xmisc/R/xdata.R
## 
## ************************************************************************

get_DEFAULT__pcThreads <- function(){
  ## 'mc.cores' > 1 is not supported on Windows
  getOption("mc.cores", 1L)
}

is.data.table.interval <- function(
  x
  ){
  "data.table.interval" %in% class(x) & is.data.table(x)
}

##' (TBA)
##'
##' Copied from package Xmisc (beta)
##' @title Convert data.table as `data.table.interval` object
##' @param x 
##' @param by 
##' @param extra.col.rm 
##' @return 
##' @author
##' @examples
##' 
dt.as_interval <- function(
  x,
  by=NULL,
  which.col.start="start",
  which.col.end="end",
  extra.col.rm=FALSE
  ){
  ret <- copy(x)
  
  if (is.data.table.interval(ret)){
    return(ret)
  }
  
  if (missing(by)){
    by <- NULL
  }
  
  if (!is.data.table(ret)){
    setDT(ret)
  }

  setnames(ret,which.col.start,"start")
  setnames(ret,which.col.end,"end")
  
  if (extra.col.rm){
    ret <- ret[,c(by,"start","end"),with=FALSE]    
  } else {
    ret <- ret[,c(by,"start","end",colnames(ret)[!colnames(ret) %in% c(by,"start","end")]),with=FALSE]
  }
  ret[,start:=as.numeric(start)]
  ret[,end:=as.numeric(end)]
  class(ret) <- unique(c("data.table.interval",class(ret)))
  return(ret)
}


  

##' (TBA)
##'
##' 
##' @title Sort data.table.interval
##' @param x 
##' @param by 
##' @return data.table
##' @author Xiaobei Zhao
dt.sort_interval <- function(
  x,
  by=NULL,
  ...
  ){
  if (missing(by)){
    by <- NULL
  }
  x2 <- dt.as_interval(x,by,...)
  .cols <- c(by,"start","end")
  .order <- c({if (is.null(by)) NULL else rep_len(1,length(by))},1,1)
  setorderv(x2,.cols,.order)
  x2
}



##' (TBA)
##'
##' 
##' @title Merge data.table.interval
##' @param x 
##' @param y 
##' @param by character to grouped by, i.e. "refid", "chr"
##' @return data.table
##' @author Xiaobei Zhao
##' @examples
##' 
dt.merge_interval <- function(x,y,by=NULL,...){

  ##
  x2 <- dt.sort_interval(x,by,...)
  y2 <- dt.sort_interval(y,by,...)
  
  ##
  suppressWarnings(x2[,start1:=start+1])
  suppressWarnings(y2[,start1:=start+1])

  ##
  setkeyv(x2,c(by,'start1','end'))
  setkeyv(y2,c(by,'start1','end'))
  ans <- foverlaps(x2, y2, type="any")
  setnames(ans,c('start','end','start1'),paste0(c('start','end','start1'),'.y'))
  setnames(ans,c('i.start','i.end','i.start1'),paste0(c('start','end','start1'),'.x'))

  ##
  ans <- ans[, `:=`(start=pmax(start1.x, start1.y)-1, end=pmin(end.x, end.y))]
  ans <- ans[, `:=`(start1.x=NULL,  start1.y=NULL)][start < end] ## start1 <= end
  ans <- dt.sort_interval(ans,by,extra.col.rm=FALSE)
    
  ##
  tmp.names <- c("start.x","end.x","start.y","end.y")
  ans <- ans[,c(colnames(ans)[!colnames(ans)%in%tmp.names],tmp.names),with=FALSE]
  ans[,size:=end-start]
  return(ans)
}



.dt.score_interval <- function(
  data.region,
  data.score,

  ##
  by,
  which.col.start,
  which.col.end,
  which.col.score,
  which.col.sampleName,
  which.col.regionName,

  ##
  FUNC.SCORE,
  score.weighted,

  ##
  by.coord,
  
  ##
  pcThreads
  ){


  
  .colnames.region <- c(by,which.col.start,which.col.end,which.col.regionName)
  if (!all(.colnames.region %in% colnames(data.region))){
    sapply(.colnames.region,stop.not.in.vector,y=colnames(data.region))
  }

  .colnames.score <- c(by,which.col.start,which.col.end,which.col.score,which.col.sampleName)
  if (!all(.colnames.score %in% colnames(data.score))){
    sapply(.colnames.score,stop.not.in.vector,y=colnames(data.score))
  }

  ## -- ##
  ## printme(list(
  ##   data.region=head(data.region),
  ##   data.score=head(data.score),
  ##   which.col.score=which.col.score,
  ##   which.col.sampleName=which.col.sampleName,
  ##   which.col.regionName=which.col.regionName
  ##   ), ".dt.score_interval")
  
  ## 
  which.col.size <- "size"

  ##
  if (is.null(data.region)){
    data.region <- data.table()
  } 
  
  if (is.null(data.score)){
    data.score <- data.table()
  }
  
  ##
  if (!nrow(data.score)){ ## empty `data.score`
    ret <- copy(data.region)
    if (nrow(data.region)){## empty `data.region`
      ret <- cbind(ret,data.table(..score..=character()))
    } else {
      ret[["..score.."]] <- NA_character_
    }
  } else {      
    tmp <- dt.merge_interval(data.region,data.score,by=by,extra.col.rm=FALSE)
    ## ## printme(tmp,"dt.score_interval | .dt.score_interval")
    ## ## printme(str(tmp),"dt.score_interval | .dt.score_interval")
    ## stop(save(file="~/save.RData",list=ls(all.names=TRUE))) ## ## load(file="~/save.RData")

    if (by.coord){
      .by.region <- c(which.col.regionName,by,"start.x","end.x")
    } else {
      .by.region <- c(which.col.regionName)      
    }
    
    ##
    if (!score.weighted){
      .sd <- tmp[,c(.SD,.SD[,list(..score..=as.character(FUNC.SCORE(.SD[[which.col.score]])))]),by=c(which.col.sampleName,.by.region)]
    } else {
      .sd <- tmp[,c(.SD,.SD[,list(..score..=as.character(FUNC.SCORE(.SD[[which.col.score]],.SD[[which.col.size]])))]),by=c(which.col.sampleName,.by.region)]      
    }
    
    ret <- .sd[!duplicated(.sd[,c(which.col.sampleName,.by.region),with=FALSE]),c(which.col.sampleName,.by.region,"..score.."),with=FALSE]
    if (by.coord){
      setnames(ret,c("start.x","end.x"),c("start","end"))
    }

    ##
    ret <- ret[!is.na(ret[[which.col.regionName]])]
    
  }

  ## -- ##

  ## ## printme(str(ret),"dt.score_interval | `.dt.score_interval`")
  if (by.coord){
    setnames(ret,"start",which.col.start)
    setnames(ret,"end",which.col.end)
  }

  
  ## 
  ret[["..score.."]] <- type.convert(ret[["..score.."]],as.is=TRUE)
  setnames(ret,"..score..",which.col.score)
  ## ## printme(str(ret),"dt.score_interval | `.dt.score_interval`")  
  
  
  return(ret)
}


FUNC__by <- function(
  e,
  by
  ){
  as.character(e[[by[1]]])
}


FUNC__by.class <- function(
  e,
  by
  ){
  class(e[[by[1]]])
}


dl.get_by.comm <- function(
  x.list,
  by=c("chr")
  ){
  .by.list <- lapply(x.list,FUNC__by,by=by)
  .by.comm <- Reduce(intersect,.by.list)  
  ret <- .by.comm
  return(ret)
}


dl.get_by.is.character <- function(
  x.list,
  by=c("chr")
  ){
  .by.is.character <- lapply(x.list,function(e){is.character(e[[by[1]]])})
  ret <- unlist(.by.is.character)
  return(ret)
}


dl.get_by.is.factor <- function(
  x.list,
  by=c("chr")
  ){
  .by.is.factor <- lapply(x.list,function(e){is.factor(e[[by[1]]])})
  ret <- unlist(.by.is.factor)
  return(ret)
}



dl.common_by <- function(
  x.list,
  by=c("chr")  
  ){
  .by.comm <- dl.get_by.comm(x.list=x.list,by=by)
  ret <- lapply(x.list,function(e){e[FUNC__by(e,by=by) %in% .by.comm]})
  return(ret)  
}



dl.factor_by <- function(
  x.list,
  by=c("chr")  
  ){
  .by.comm <- dl.get_by.comm(x.list=x.list,by=by)
  .by.is.factor <- dl.get_by.is.factor(x.list=x.list,by=by)
  ## ## .by.is.character <- dl.get_by.is.character(x.list=x.list,by=by)

  if (!all(.by.is.factor)){
    .levels <- gtools::mixedsort(.by.comm)
    ret <- lapply(x.list,function(e){e[[by[1]]] <- factor(as.character(e[[by[1]]]),levels=.levels); return(e)})
  } else {
    ret <- x.list
  }
  return(ret)
}
  


##' To score `data.region` using `data.score` by subgrouped regions.
##'
##' In computational biology, it scores a set of regions across chromosomes using a set of scores by overlapping genomic coordinates.
##' @title To score `data.region` using `data.score`
##' @param data.region a data table of regions to score. i.e. a data table of region information by feature (e.g. gene).
##' @param data.score a data table of score information over regions. i.e. a data table of score information over regions by sample.
##' @param by the column, within each level of which `data.region` and `data.score` need be processed individually. (e.g. chromosome)
##' @param which.col.start the column for the star coordinate of a region for overlapping. (0-based)
##' @param which.col.end the column for the end coordinate of a region for overlapping. (1-based)
##' @param which.col.score the column for score
##' @param which.col.sampleName the column for sample identifiers
##' @param which.col.regionName the column for region (feature, e.g. gene) identifiers
##' @param FUNC.SCORE a function to compute a single score given a vector of scores. (e.g. mean or median)
##'        see \code{\link{get_FUNC.SCORE_mean}}, \code{\link{get_FUNC.SCORE_median}}
##' @param pcThreads the processing threads. (pcThreads > 1 is not supported on Windows)
##' @return a data table with region name, coordinates, sample name and scores over a region per sample.
##' @author Xiaobei Zhao
##' @examples
dt.score_interval <- function(
  data.region,
  data.score,
  ##
  by=c("chr"),
  which.col.start="start",
  which.col.end="end",
  which.col.score="score",
  which.col.sampleName="sampleName",
  which.col.regionName="regionName",

  ##
  FUNC.SCORE=NULL,
  score.weighted=FALSE,

  ##
  by.coord=TRUE,
  
  ##
  pcThreads=get_DEFAULT__pcThreads()
  ){
  
  ## printme(list(
  ##   data.region=head(data.region),
  ##   data.score=head(data.score),
  ##   which.col.score=which.col.score,
  ##   which.col.sampleName=which.col.sampleName,
  ##   which.col.regionName=which.col.regionName
  ##   ), "dt.score_interval")
  

  ## 
  options(scipen=9)
  
  if (is.null(FUNC.SCORE)){
    FUNC.SCORE <- get_FUNC.SCORE_vconcat
  }
  
  ## ## stop(save(file="~/save.RData",list=ls(all.names=TRUE))) # load(file="~/save.RData")
  
  
  
  ## if (is.null(data.region)){ ##TBA
  ##   data.region <-  unoverlap_data.region(data.score[,c(by,which.col.start,which.col.end),with=FALSE],by=by,pcThreads=pcThreads)
  ## }
  
  ## -- ##
  data.list0 <- data.list <- list(data.region,data.score)
  data.list <- dl.common_by(x.list=data.list,by=by)
    
  data.region <- data.list[[1]]
  data.score <- data.list[[2]]

  ## -- ##
  data.score <- prep_data.score(
    data.score,
    which.col.chr=by, 
    which.col.start=which.col.start,
    which.col.end=which.col.end,
    which.col.sampleName=which.col.sampleName,
    which.col.score=which.col.score
    )
  which.col.sampleName <- attr(data.score,"which.col.sampleName")
  which.col.score <- attr(data.score,"which.col.score")

  data.region <- prep_data.region(
    data.region,
    which.col.chr=by, 
    which.col.start=which.col.start,
    which.col.end=which.col.end,
    which.col.regionName=which.col.regionName
    )
  which.col.regionName <- attr(data.region,"which.col.regionName")

  
  
  ##
  .by <- .which.col.chr <- c("chr")
  .which.col.start <- "start"
  .which.col.end <- "end"
  .which.col.score <- "score"
  .which.col.sampleName <- "sampleName"
  .which.col.regionName <- "regionName"
  
  data.region <- dt.sort_interval(data.region,by=.which.col.chr,which.col.start=.which.col.start,which.col.end=.which.col.end)
  data.score <- dt.sort_interval(data.score,by=.which.col.chr,which.col.start=.which.col.start,which.col.end=.which.col.end)

  
  ## -- ##
  data.list <- list(data.region,data.score)
  .by.is.character <- dl.get_by.is.character(x.list=data.list,by=.by)
  data.list <- dl.factor_by(x.list=data.list,by=.by)
  data.region <- data.list[[1]]
  data.score <- data.list[[2]]
  
  
  ## -- ##
  if ("..score.." %in% colnames(data.region) | "..score.." %in% colnames(data.score)){
    message(str(data.region))
    message(str(data.score))
    stop('`..score..` is a reserved column name.')
  }

  
  if (pcThreads==1) {
    ret <- .dt.score_interval(
      data.region=data.region,
      data.score=data.score,
      by=.by,
      which.col.start=.which.col.start,
      which.col.end=.which.col.end,
      which.col.score=.which.col.score,
      which.col.sampleName=.which.col.sampleName,
      which.col.regionName=.which.col.regionName,
      FUNC.SCORE=FUNC.SCORE,
      score.weighted=score.weighted,
      by.coord=by.coord,
      pcThreads=pcThreads
      )
  } else {
    data.region.list <- split(data.region,as.list(data.region[,.by[1],with=FALSE]))##XB
    data.score.list <- split(data.score,as.list(data.score[,.by[1],with=FALSE]))##XB

    ##  
    tmp <- parallel::mclapply(
      names(data.region.list),
      function(a){
        data.region <- data.region.list[[a]]
        data.score <- data.score.list[[a]]
        message("dt.score_interval | ",a);
        .dt.score_interval(
          data.region=data.region,
          data.score=data.score,
          by=.by,
          which.col.start=.which.col.start,
          which.col.end=.which.col.end,
          which.col.score=.which.col.score,
          which.col.sampleName=.which.col.sampleName,
          which.col.regionName=.which.col.regionName,
          FUNC.SCORE=FUNC.SCORE,
          score.weighted=score.weighted,
          by.coord=by.coord,
          pcThreads=1
          )
      },
      mc.cores=pcThreads
      )

    ret <- do.call(rbind,tmp)
  }
  
  if (by.coord){
    if (any(.by.is.character)){
      ret[[.by[1]]] <- as.character(ret[[.by[1]]])
    }
  }

  
  setnames(ret,.which.col.regionName,which.col.regionName)
  setnames(ret,.which.col.sampleName,which.col.sampleName)
  setnames(ret,.which.col.score,which.col.score)
  
  ## printme(str(ret),"dt.score_interval")
  return(ret)
}


## ------------------------------------------------------------------------
## 
## ------------------------------------------------------------------------



dt.mscore_interval <- function(
  data.region,
  data.score,
  ##
  by=c("chr"),
  by.x=by,
  by.y=by,
  ##
  which.col.start="start",
  which.col.end="end",
  which.col.score="score",
  which.col.sampleName="sampleName",
  which.col.regionName="regionName",

  ##
  FUNC.SCORE=NULL,

  ##
  by.coord=TRUE,
  
  ##
  pcThreads=get_DEFAULT__pcThreads()
  ){
  
  .data.region <- copy(data.region)
  .data.score <- copy(data.score)
  
  setnames(.data.region,by.x,by)
  setnames(.data.score,by.y,by)

  by <- by[!by %in% c(which.col.start,which.col.end,which.col.score,which.col.sampleName,which.col.regionName)]

  ## ## printme(list(data.region=head(data.region),by.x=by.x,data.score=head(data.score),by.y=by.y),"dt.mscore_interval")
  ##
  ## ## stop(save(file="~/.save.RData",list=ls(all.names=TRUE))) # load(file="~/.save.RData")
  
  ## 
  n <- 0
  for (e in which.col.score){
    n <- n+1
    itmp <- dt.score_interval(
      data.region=.data.region,
      data.score=.data.score,

      ##
      by=by,
      which.col.start=which.col.start,
      which.col.end=which.col.end,
      which.col.score=e,
      which.col.sampleName=which.col.sampleName,
      which.col.regionName=which.col.regionName,
      
      ##
      FUNC.SCORE=FUNC.SCORE,
      
      ##
      by.coord=by.coord,
      
      ##
      pcThreads=pcThreads
      )
    if (n==1){
      tmp <- copy(itmp)
    } else {
      tmp <- cbind(tmp,itmp[,ncol(itmp),with=FALSE])
    }
  }
  ret <- tmp
  ## printme(str(ret),"dt.mscore_interval")
  return(ret)
}





## ************************************************************************
## 
## ************************************************************************
