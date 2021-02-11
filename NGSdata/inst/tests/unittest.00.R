
## ************************************************************************
## 
## 
## 
## (c) Xiaobei Zhao
## 
## Thu Feb 11 04:59:23 CST 2021 -0600 (Week 06)
## 
## 
## Reference: 
## 
## To Run:
## Rscript -e "require(Xmisc); require(NGSdata); source(\"unittest.00.R\"); unittest.00(); "
## 
## ************************************************************************


unittest.00 <- function(
  ){
  unittest.00__UnitTest()
}


unittest.00__UnitTest <-
  function()
{
  
  pkg <- 'NGSdata'
  test.obj <- UnitTest$new(pkg=pkg)
  test.obj$defineme()
  logme(test.obj)
}

