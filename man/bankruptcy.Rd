\name{bankruptcy}
\alias{bankruptcy}
\docType{data}

\title{Bankruptcy data}

\description{ The data set contain the ratio of retained earnings (RE) to total assets, and the ratio of earnings before interests and taxes (EBIT) to total assets of 66 American firms recorded in the form of ratios. Half of the selected firms had filed for bankruptcy. }

\usage{data(bankruptcy)}

\format{A data frame with the following variables:

\describe{
	\item{Y}{The status of the firm: \code{0} bankruptcy or \code{1} financially sound;}
	\item{RE}{Ratio of retained earnings to total assets;}
	\item{EBIT}{Ratio of earnings before interests and taxes to total assets}

}}
\references{
Altman E.I. (1968) Financial ratios, discriminant analysis and the prediction of corporate bankruptcy, \emph{J Finance} \bold{23}(4): 589-609
}

\examples{

data(bankruptcy)

}

\keyword{data sets}
