#' Swedish municipalities data
#'
#' The panel data set consists of 265 Swedish municipalities and
#' covers 9 years (1979-1987).
#'
#' @format The variables are:
#' \describe{
#'   \item{id}{ID number for municipality}
#'   \item{year}{Year}
#'   \item{expenditures}{Total expenditures}
#'   \item{revenues}{Total own-source revenues}
#'   \item{grants}{Intergovernmental grants received by the municipality}
#' }
#'
#' Total expenditures contains both capital and current expenditures.
#'
#' Expenditures, revenues, and grants are expressed in million SEK. The
#' series are deflated and in per capita form. The implicit deflator is a
#' municipality-specific price index obtained by dividing total local
#' consumption expenditures at current prices by total local consumption
#' expenditures at fixed (1985) prices.
#'
#' The data are gathered by Statistics Sweden and obtained from
#' Financial Accounts for the Municipalities (Kommunernas Finanser).
#' @source \url{http://qed.econ.queensu.ca/jae/2000-v15.4/dahlberg-johansson/}
#' @references M. Dahlberg and E. Johansson (2000) "An examination of the dynamic behavior of local governments using GMM bootstrapping methods", \emph{Journal of Applied Econometrics}, \strong{15}(4), 401-416, \url{http://www.jstor.org/stable/2678589}.
"Dahlberg"

#' Employment UK data
#'
#' This data set contains labor demand data from a panel of firms in the United Kingdom. The panel is unlanced.
#'
#' @format The variables are:
#' \describe{
#'   \item{c1}{Record ID}
#'   \item{ind}{Firm index}
#'   \item{year}{Year}
#'   \item{emp}{Employment}
#'   \item{wage}{Wage}
#'   \item{cap}{Capital}
#'   \item{indoutpt}{Industrial output}
#'   \item{n, w, k, ys}{Logs of variables}
#'   \item{rec}{Record number}
#'   \item{yearm1}{Lagged year}
#'   \item{id}{ID}
#'   \item{nL1, nL2, wL1, kL1, kL2, ysL1, ysL2}{Lags of log variables}
#'   \item{yr1976 - yr1984}{Time dummies}
#' }
#'
#' @source \url{http://www.stata-press.com/data/r13/abdata.dta}
#' @references Arellano, M. and Bond, S. (1991) "Some tests of specification for panel data: Monte Carlo evidence and an application to employment equations", \emph{The Review of Economic Studies}, \strong{58}(2), 227-297, \href{http://dx.doi.org/10.2307/2297968}{doi: 10.2307/2297968}
"abdata"

#' Cigar data
#'
#' This panel data set consists of 46 U.S. States over the period 1963-1992.
#'
#' @format The variables are:
#' \describe{
#'   \item{state}{State abbreviation}
#'   \item{year}{Year}
#'   \item{price}{Price per pack of cigarettes}
#'   \item{pop}{Population}
#'   \item{pop16}{Population above the age of 16.}
#'   \item{cpi}{Consumer price index with (1983=100}
#'   \item{ndi}{Per capita disposable income}
#'   \item{sales}{Cigarette sales in packs per capita}
#'   \item{pimin}{Minimum price in adjoining states per pack of cigarettes}
#'
#' }
#'
#' All variables all also available as logs.
#'
#' @source \url{http://www.wiley.com/legacy/wileychi/baltagi/supp/Cigar.txt}
#' @references
#' Baltagi, B.H. and D. Levin (1992) "Cigarette taxation: raising revenues and reducing consumption", \emph{Structural Change and Economic Dynamics}, \strong{3}(2), 321-335, \href{http://dx.doi.org/10.1016/0954-349X(92)90010-4}{doi:10.1016/0954-349X(92)90010-4}.
#'
#' Baltagi, B.H., J.M. Griffin and W. Xiong (2000) "To pool or not to pool: homogeneous versus heterogeneous estimators applied to cigarette demand", \emph{Review of Economics and Statistics}, \strong{82}(1), 117-126, \href{http://dx.doi.org/10.1162/003465300558551}{doi:10.1162/003465300558551}.
#'
#' Baltagi, B.H. (2013) "Econometric analysis of panel data", 5th edition, \emph{John Wiley and Sons}, \url{http://www.wiley.com/WileyCDA/WileyTitle/productCd-EHEP003191.html}

#' Cigar
"Cigar"

#' NLS Work 2 data
"nlswork2"

#' Dahlberg results example 1
"ex1_dahlberg_data"

#' Dahlberg bootstrap results example 1
"ex1_dahlberg_data_bs"

#' NLS Work 2 bootstrap results example 2
"ex2_nlswork2_data_bs"

#' Example results for Employment UK data
"ex3_abdata"