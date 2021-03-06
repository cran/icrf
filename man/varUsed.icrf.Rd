% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varUsed.R
\name{varUsed.icrf}
\alias{varUsed.icrf}
\alias{varUsed}
\title{'Variables used in' an \code{icrf} ensemble}
\usage{
varUsed(x, ...)

\method{varUsed}{icrf}(x, by.tree = FALSE, count = TRUE, ...)
}
\arguments{
\item{x}{'an object of class' \code{icrf}.}

\item{...}{not used.}

\item{by.tree}{'Should the list of variables used be broken down by trees in the forest?'}

\item{count}{'Should the frequencies that variables appear in trees be returned?'}
}
\value{
'A vector containing number of nodes for the trees' in the icrf object.

'If \code{count=TRUE} and \code{by.tree=FALSE}, an integer vector containing
frequencies that variables are used in the forest. If \code{by.tree=TRUE}, a matrix is returned,
breaking down the counts by tree (each column corresponding to one tree and each row to a variable).'

'If \code{count=FALSE} and \code{by.tree=TRUE}, a list of integer indices is returned giving the
variables used in the trees, else if \code{by.tree=FALSE}, a vector of integer indices giving
the variables used in the entire forest.'
}
\description{
'Find out which predictor variables are actually used in' the returned forest of
the \code{icrf}. The returned forest depends on the \code{returnBest} argument of the
\code{icrf} function; It is either the last forest, when \code{returnBest = FALSE} or
the the best forest, when \code{returnBest = TRUE}. (Quoted statements are from
\code{randomForest} by Liaw and Wiener unless otherwise mentioned.)
}
\examples{
# rats data example.
# Note that this is a toy example. Use a larger ntree and nfold in practice.
data(rat2)
\donttest{
 set.seed(1)
 rats.icrf <-
   icrf(~ dose.lvl + weight + male + cage.no, data = rat2,
        data.type = "currentstatus", currentstatus.label = c("survtime", "tumor"),
        returnBest = TRUE, ntree=10, nfold=3)
 varUsed(rats.icrf)
}
\dontshow{
 set.seed(1)
 rats.icrf <-
   icrf(~ dose.lvl + weight + male + cage.no, data = rat2,
        data.type = "currentstatus", currentstatus.label = c("survtime", "tumor"),
        returnBest = TRUE, ntree=2, nfold=2)
 varUsed(rats.icrf)
}

}
\references{
\href{https://arxiv.org/abs/1912.09983}{Cho H., Jewell N. J., and Kosorok M. R. (2020+). "Interval censored
 recursive forests"}
}
\author{
Hunyong Cho, Nicholas P. Jewell, and Michael R. Kosorok.
}
