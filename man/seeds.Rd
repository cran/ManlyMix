\name{seeds}
\alias{seeds}
\docType{data}
\encoding{UTF-8}
\title{Wheat kernel Data}
\description{The examined group comprised kernels belonging to three different varieties of wheat: Kama, Rosa and Canadian, 70 elements each, randomly selected for the experiment. High quality visualization of the internal kernel structure was detected using a soft X-ray technique. Studies were conducted using combine harvested wheat grain originating from experimental fields, explored at the Institute of Agrophysics of the Polish Academy of Sciences in Lublin. 
}
\usage{data(seeds)}
\format{ 
  A data frame with 210 observations on the following 7 variables.
  \describe{
   \item{V1}{Area A;}
   \item{V2}{Perimeter P;}
   \item{V3}{Compactness;}
   \item{V4}{Length of kernel;}
   \item{V5}{Width of kernel;}
   \item{V6}{Asymmetry coefficient;}
   \item{V7}{Length of kernel groove;}
   \item{V8}{Seed species: \code{1}, \code{2}, \code{3}}
   }
}
\references{
M. Charytanowicz, J. Niewczas, P. Kulczycki, P.A. Kowalski, S. Lukasik, S. Zak (2010), 
\emph{A Complete Gradient Clustering Algorithm for Features Analysis of X-ray Images}.
Information Technologies in Biomedicine, Ewa Pietka, Jacek Kawa, Springer-Verlag, Berlin-Heidelberg.
}

\examples{

data(seeds)

}
\keyword{datasets}
