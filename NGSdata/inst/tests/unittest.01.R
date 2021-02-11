
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Thu Feb 11 05:11:46 CST 2021 -0600 (Week 06)
## 
## 
## Reference: 
## Rscript -e "require(Xmisc); require(NGSdata); source(\"unittest.01.R\"); unittest.01(); "
## 
## ************************************************************************



unittest.01 <- function(){
  options(width=150)
  ## unittest.01__test()
  ## unittest.02__test()
  unittest.02__example()
} 



unittest.01__test <- function(){
  printme(get_toy__data.score())
  printme(get_toy__data.region())  
}

unittest.02__test <- function(){
  
  ##
  .data.score <- read.data(get_toy__data.score())
  .data.region <- read.data(get_toy__data.region())

  ##
  .data.score <- prep_data.score(
    .data.score,
    which.col.sampleName="sampleID",
    which.col.score="sgmt.value"
    )

  .data.region <- prep_data.region(
    .data.region,
    which.col.regionName="featureID"
    )

  ##
  .data.score <- dt.as_interval(.data.score,by="chr")
  .data.region <- dt.as_interval(.data.region,by="chr")
  
  ## ##
  ## .data.score <- dt.sort_interval(.data.score,by=NULL)
  ## .data.region <- dt.sort_interval(.data.region,by=NULL)
  
  ##
  .data.score <- dt.sort_interval(.data.score,by="chr")
  .data.region <- dt.sort_interval(.data.region,by="chr")

  ##
  tmp.merge <- dt.merge_interval(
    .data.region,
    .data.score,
    by="chr"
    )
  
  tmp.score <- dt.score_interval(
    .data.region,
    .data.score,
    by="chr",
    FUNC.SCORE=get_FUNC.SCORE_mean,
    pcThreads=1
    )

  tmp.mat <- data.as.matrix(
    tmp.score
    )
  
  ##
  printme((.data.score))
  printme((.data.region))
  printme(tmp.merge)
  printme(tmp.score)
  printme(tmp.mat)

  ##
  printme(packageVersion("Xmisc"))
  printme(packageVersion("NGSdata"))
  
}




unittest.02__example <- function(){

  data.score <- read.data(get_toy__data.score())
  data.region <- read.data(get_toy__data.region())
  
  ##
  .data.score <- prep_data.score(
    data.score,
    which.col.chr="chr",
    which.col.start="start",
    which.col.end="end",
    which.col.sampleName="sampleID",
    which.col.score="sgmt.value"
    )

  .data.region <- prep_data.region(
    data.region,
    which.col.chr="chr",
    which.col.start="start",
    which.col.end="end",
    which.col.regionName="featureID"
    )

  ##
  tmp.score <- dt.score_interval(
    data.region=.data.region,
    data.score=.data.score,
    by="chr",
    FUNC.SCORE=function(x){mean(x,na.rm=TRUE)},
    pcThreads=2
  )
  
  tmp.mat <- data.as.matrix(
    tmp.score
    )
  
  printme((.data.score))
  printme((.data.region))
  printme(tmp.mat)
  invisible()
}
