% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.icrf.R
\name{predict.icrf}
\alias{predict.icrf}
\title{icrf predictions}
\usage{
\method{predict}{icrf}(
  object,
  newdata,
  predict.all = FALSE,
  proximity = FALSE,
  nodes = FALSE,
  smooth = TRUE,
  ...
)
}
\arguments{
\item{object}{an object of \code{icrf} class generated by the function \code{icrf}.}

\item{newdata}{'a data frame or matrix containing new data. (Note: If not given,
the predicted survival estimate of the training data set in the \code{object} is returned.)'}

\item{predict.all}{'Should the predictions of all trees be kept?'}

\item{proximity}{'Should proximity measures be computed?'}

\item{nodes}{'Should the terminal node indicators (an n by ntree matrix)
be returned? If so, it is in the "nodes" attribute of the returned object.'}

\item{smooth}{Should smoothed curve be returned?}

\item{...}{'not used currently.'}
}
\value{
A matrix of predicted survival probabilities is returned where the rows represent
the observations and the columns represent the time points.
'If predict.all=TRUE, then the returned object is a list of two components:
\code{aggregate}, which is the vector of predicted values by the forest,
and \code{individual}, which is a matrix where each column contains prediction
by a tree in the forest.' The forest is either the last forest or the best forest
as specified by \code{returnBest} argument in \code{icrf} function.


'If \code{proximity=TRUE}, the returned object is a list with two components:
\code{pred} is the prediction (as described above) and
\code{proximity} is the proximitry matrix.'


'If \code{nodes=TRUE}, the returned object has a "nodes" attribute,
which is an \code{n} by \code{ntree} matrix, each column containing the
node number that the cases fall in for that tree.'
}
\description{
Prediction method of test data using interval censored recursive forest.
(Quoted statements are from
\code{randomForest} by Liaw and Wiener unless otherwise mentioned.)
}
\examples{
# rats data example
# Note that this is a toy example. Use a larger ntree and nfold in practice.
library(survival)  # for Surv()
data(rat2)
set.seed(1)
samp <- sample(1:dim(rat2)[1], 200)
rats.train <- rat2[samp, ]
rats.test <- rat2[-samp, ]
L = ifelse(rats.train$tumor, 0, rats.train$survtime)
R = ifelse(rats.train$tumor, rats.train$survtime, Inf)
\donttest{
 set.seed(2)
 rats.icrf.small <-
   icrf(survival::Surv(L, R, type = "interval2") ~ dose.lvl + weight + male + cage.no,
        data = rats.train, ntree = 10, nfold = 3, proximity = TRUE)

 # predicted survival curve for the training data
 predict(rats.icrf.small)
 predict(rats.icrf.small, smooth = FALSE) # non-smoothed

 # predicted survival curve for new data
 predict(rats.icrf.small, newdata = rats.test)
 predict(rats.icrf.small, newdata = rats.test, proximity = TRUE)

 # time can be extracted using attr()
 newpred = predict(rats.icrf.small, newdata = rats.test)
 attr(newpred, "time")

 newpred2 = predict(rats.icrf.small, newdata = rats.test, proximity = TRUE)
 attr(newpred2$predicted, "time")
}
\dontshow{
 set.seed(2)
 rats.icrf.small <-
   icrf(survival::Surv(L, R, type = "interval2") ~ dose.lvl + weight + male + cage.no,
        data = rats.train, ntree = 2, nfold = 2, proximity = TRUE)

 # predicted survival curve for the training data
 predict(rats.icrf.small)
 predict(rats.icrf.small, smooth = FALSE) # non-smoothed

 # predicted survival curve for new data
 predict(rats.icrf.small, newdata = rats.test)
 predict(rats.icrf.small, newdata = rats.test, proximity = TRUE)

 # time can be extracted using attr()
 newpred = predict(rats.icrf.small, newdata = rats.test)
 attr(newpred, "time")

 newpred2 = predict(rats.icrf.small, newdata = rats.test, proximity = TRUE)
 attr(newpred2$predicted, "time")
}


}
\references{
\href{https://arxiv.org/abs/1912.09983}{Cho H., Jewell N. J., and Kosorok M. R. (2020+). "Interval censored
 recursive forests"}
}
\author{
Hunyong Cho, Nicholas P. Jewell, and Michael R. Kosorok.
}
