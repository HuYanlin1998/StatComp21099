\name{StatComp21099-package}
\alias{StatComp21099-package}
\alias{StatComp21099}
\docType{package}
\title{
  Test difference between two survival curves.
}
\description{
  Test if there are difference between two survival curves,using Fleming-Harrington and Peto&Peto method,based on both big sample approximation and permutation test.
}
\details{
  This package includes two major functions :"surtest"and "gibbsC",the first is to test difference between two survival curves,and the second is to generate random samples using gibbs sampling ,which is a Cfunction.More details are in comments.
}
\author{
Yanlin Hu

Maintainer: Yanlin Hu <hyl11@mail.ustc.edu.cn>
}

\keyword{ package }

\examples{
  \dontrun{
     ## tableCaculation(aml1)
     ## SurvFuncEsti(aml1)
     ## survtest(aml1,999,"FH",1,2)
     ## survtest(aml1,999,"PP",modified=FALSE)
  }
}
