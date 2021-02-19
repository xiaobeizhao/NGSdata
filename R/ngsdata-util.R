
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


##' .. content for \description{} (no empty lines) ..
##'
##' Copied from package Xmisc (beta)
##' @title 
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
  if (missing(by)){
    by <- NULL
  }
  ret <- copy(x)
  
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
  return(ret)
}


  

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title dt.sort_interval
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



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title dt.merge_interval
##' @param x 
##' @param y 
##' @param by character to grouped by, i.e. "chr", "refid"
##' @return data.table
##' @author Xiaobei Zhao
##' @examples
##' 
dt.merge_interval <- function(x,y,by=NULL,...){

  ##
  x2 <- dt.sort_interval(x,by,...)
  y2 <- dt.sort_interval(y,by,...)
  
  ##
  x2[,start1:=start+1]
  y2[,start1:=start+1]

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
  return(ans)
}





##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title To score `data.region` using `data.score`
##' @param data.region 
##' @param data.score 
##' @param by 
##' @param which.col.start 
##' @param which.col.end 
##' @param which.col.score 
##' @param FUNC.SCORE 
##' @param pcThreads 
##' @param keepname 
##' @param as.dt 
##' @return 
##' @author 
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

  ##
  pcThreads
  ){
  
  
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

    .ddply <- plyr::ddply(
      tmp,
      c(which.col.sampleName,which.col.regionName,by,"start.x","end.x"),
      function(df){
        ## printme(df)
        ## printme(FUNC.SCORE)
        .score <- df[[which.col.score]]
        .score <- FUNC.SCORE(.score)
        .score <- as.character(.score)
        .df <- data.table(..score..=.score)
        .df
      }
      )
    
    setDT(.ddply)
    if (!nrow(.ddply)){
      .ddply <- cbind(.ddply[,c(by,"start.x","end.x"),with=FALSE],data.table("..score.."=character(),stringsAsFactors=FALSE))
    }
    ## printme(colnames(.ddply))
    ## printme(str(.ddply))

    .merge <- merge(
      data.region[,unique(c(c(which.col.regionName,by,"start","end"),colnames(data.region)[!colnames(data.region) %in% c(which.col.score,colnames(.ddply))])),with=FALSE],
      .ddply,
      by.x=c(which.col.regionName,by,"start","end"),by.y=c(which.col.regionName,by,"start.x","end.x"),
      all=TRUE,
      allow.cartesian=TRUE,
      sort=FALSE
      )
    ##
    ret <- .merge[!is.na(.merge[[which.col.regionName]])]
  }

  ## -- ##
  setDT(ret)
  setnames(ret,"start",which.col.start)
  setnames(ret,"end",which.col.end)

  ## 
  setnames(ret,"..score..",which.col.score)
  ## printme(str(ret),"dt.score_interval | `.dt.score_interval`")
  ## printme(str(ret[[which.col.score]]))
  
  ret[[which.col.score]] <- type.convert(ret[[which.col.score]],as.is=TRUE)
  setDF(ret)
  
  return(ret)
}








##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param data.region 
##' @param data.score 
##' @param by 
##' @param which.col.score 
##' @param FUNC.SCORE 
##' @return 
##' @author Xiaobei Zhao
##' @examples
##' 
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

  pcThreads=get_DEFAULT__pcThreads()
  ){
  
  ## 
  options(scipen=9)
  
  if (is.null(FUNC.SCORE)){
    FUNC.SCORE <- get_FUNC.SCORE_vconcat
  }
  
  
  ## if (is.null(data.region)){ ##TBA
  ##   data.region <-  unoverlap_data.region(data.score[,c(by,which.col.start,which.col.end),with=FALSE],by=by,pcThreads=pcThreads)
  ## }

  
  ## -- ##
  .by.score <- as.character(data.score[[by[1]]])
  .by.region <- as.character(data.region[[by[1]]])
  .by <- unique(c(.by.region,.by.score))
  .by <- .by[.by %in% .by.region & .by %in% .by.score]
  
  data.region <- data.region[as.character(data.region[[by[1]]]) %in% .by]
  data.score <- data.score[as.character(data.score[[by[1]]]) %in% .by]
  
  data.region <- dt.sort_interval(data.region,by,which.col.start=which.col.start,which.col.end=which.col.end)
  data.score <- dt.sort_interval(data.score,by,which.col.start=which.col.start,which.col.end=which.col.end)

  
  ## -- ##
  
  by1.class <- class(data.region[[by[1]]])
  if (!is.factor(data.score[[by[1]]]) | !is.factor(data.region[[by[1]]])){
    .levels <- gtools::mixedsort(.by)
    data.region[[by[1]]] <- factor(data.region[[by[1]]],levels=.levels)
    data.score[[by[1]]] <- factor(data.score[[by[1]]],levels=.levels)
  }
  
  
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
      by=by,
      which.col.start=which.col.start,
      which.col.end=which.col.end,
      which.col.score=which.col.score,
      which.col.sampleName=which.col.sampleName,
      which.col.regionName=which.col.regionName,
      FUNC.SCORE=FUNC.SCORE,
      pcThreads=pcThreads
      )
  } else {
    data.region.list <- split(data.region,as.list(data.region[,by[1],with=FALSE]))##XB
    data.score.list <- split(data.score,as.list(data.score[,by[1],with=FALSE]))##XB

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
          by=by,
          which.col.start=which.col.start,
          which.col.end=which.col.end,
          which.col.score=which.col.score,
          which.col.sampleName=which.col.sampleName,
          which.col.regionName=which.col.regionName,
          FUNC.SCORE=FUNC.SCORE,
          pcThreads=1
          )
      },
      mc.cores=pcThreads
      )

    ret <- do.call(rbind,tmp)
    ## ## .order <- do.call(order,as.list(ret[,c(by,"start","end"),with=FALSE]))
    ## ## ret <- ret[.order,]
  }
  
  if ( "character" %in% by1.class){
    ret[[by[1]]] <- as.character(ret[[by[1]]])
  }

  ## printme(str(ret),"dt.score_interval")
  return(ret)
}






## ************************************************************************
## 
## ************************************************************************
