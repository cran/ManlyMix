%  file ManlyMix/man/ais.Rd  
%  This file is a component of the package 'ManlyMix' for R
%  copyright (C) 2004-2013 Adelchi Azzalini
%---------------------
\name{ais}
\alias{ais}
\docType{data}
\encoding{UTF-8}
\title{Australian Institute of Sport data}
\description{Data on 102 male and 100 female athletes collected at the 
Australian Institute of Sport, courtesy of Richard Telford and Ross Cunningham.
}
\usage{data(ais)}
\format{ 
  A data frame with 202 observations on the following 13 variables.
  \describe{
   \item{sex}{Factor with levels: \code{female}, \code{male};}
   \item{sport}{Factor with levels: \code{B_Ball}, 
            \code{Field}, \code{Gym}, \code{Netball}, \code{Row}
            \code{Swim}, \code{T_400m}, \code{Tennis}, \code{T_Sprnt}, 
            \code{W_Polo};}
   \item{RCC}{Red cell count;}
   \item{WCC}{White cell count;}
   \item{Hc}{Hematocrit;}
   \item{Hg}{Hemoglobin;}
   \item{Fe}{Plasma ferritin concentration;}
   \item{BMI}{Body Mass Index;}
   \item{SSF}{Sum of skin folds;}
   \item{Bfat}{Body fat percentage;}
   \item{LBM}{Lean body mass;}
   \item{Ht}{Height, cm;}
   \item{Wt}{Weight, kg}
  }
}
\details{The data have been made publicly available in connection with the
book by  Cook and Weisberg (1994).}

\references{
  Cook and Weisberg (1994) \emph{An Introduction to Regression Graphics}, John Wiley & Sons, New York.
}

\examples{

data(ais)

}
\keyword{datasets}
