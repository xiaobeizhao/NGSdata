
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
## Rscript -e "require(Xmisc); require(NGSdata); source('unittest.01.R'); unittest.01(); "
## Rscript -e "require(Xmisc); require(NGSdata); source(system.file('tests','unittest.01.R',package='NGSdata',mustWork=TRUE)); unittest.01(); "
## 
## ************************************************************************



unittest.01 <- function(){
  options(width=150)
  ## unittest.01__test()
  unittest.01__example()
} 



unittest.01__test <- function(){
  printme(get_toy__data.score())
  printme(get_toy__data.region()) 
}



unittest.01__example <- function(){

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
  tmp.merge <- dt.merge_interval(
    .data.region,
    .data.score,
    by="chr"
    )
  
  ##
  tmp.score <- dt.score_interval(
    data.region=.data.region,
    data.score=.data.score,
    by="chr",
    FUNC.SCORE=function(x){mean(x,na.rm=TRUE)},
    pcThreads=get_DEFAULT__pcThreads()
  )
  
  ##
  tmp.score.weighted <- dt.score_interval(
    data.region=.data.region,
    data.score=.data.score,
    by="chr",
    FUNC.SCORE=function(x,size){sum(x*size/sum(size,na.rm=TRUE),na.rm=TRUE)},
    score.weighted=TRUE,
    pcThreads=get_DEFAULT__pcThreads()
    )

  
  ##
  tmp.mat <- data.as.matrix(
    tmp.score
    )
  
  ##
  tmp.mat.transposed <- data.as.matrix(
    tmp.score,
    regionName ~ sampleName
    )
  
  
  ##
  tmp.mat.weighted <- data.as.matrix(
    tmp.score.weighted
    )
  
  ##
  tmp.mat.weighted.transposed <- data.as.matrix(
    tmp.score.weighted,
    regionName ~ sampleName
    )
  
  ##
  printme(.data.score)
  printme(.data.region)
  printme(tmp.merge)
  printme(tmp.score)
  printme(tmp.score.weighted)
  printme(tmp.mat)
  printme(tmp.mat.transposed)
  printme(tmp.mat.weighted)
  printme(tmp.mat.weighted.transposed)
  invisible()
  
}
