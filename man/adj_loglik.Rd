% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj_loglik_methods.R
\name{adj_loglik}
\alias{adj_loglik}
\alias{coef.chantrics}
\alias{confint.chantrics}
\alias{coeftest.chantrics}
\alias{conf_region.chantrics}
\alias{conf_intervals.chantrics}
\alias{logLik.chantrics}
\alias{plot.chantrics}
\alias{print.chantrics}
\alias{summary.chantrics}
\alias{vcov.chantrics}
\title{Loglikelihood adjustments for fitted models}
\usage{
adj_loglik(x, cluster = NULL, use_vcov = TRUE, use_mle = TRUE, ...)
}
\arguments{
\item{x}{A supported fitted model object, see \strong{Supported models}}

\item{cluster}{A vector or factor indicating the cluster the corresponding
loglikelihood contribution belongs to. It is required to have the same
length as the vector returned by \code{\link[=logLik_vec]{logLik_vec()}}. If \code{cluster} is not
supplied or \code{NULL}, then it is assumed that each observation forms its own
cluster.}

\item{use_vcov}{A logical scalar. By default, the \code{\link[=vcov]{vcov()}} method for \code{x} is
used to estimate the Hessian of the independence loglikelihood, if the
function exists. Otherwise, or if \code{use_vcov = FALSE}, \code{H} is estimated
using \code{\link[stats:optim]{stats::optimHess()}} inside \code{\link[chandwich:adjust_loglik]{chandwich::adjust_loglik()}}.}

\item{use_mle}{A logical scalar. By default, the MLE from \code{x} is taken as
given, and is not reestimated. By setting \code{use_mle} to \code{FALSE}, the parameters
are reestimated in the function \code{\link[chandwich:adjust_loglik]{chandwich::adjust_loglik()}} using
\code{\link[stats:optim]{stats::optim()}}.This feature is currently for development purposes only,
may return misleading/false results and may be removed without notice.}

\item{...}{Further arguments to be passed to \code{\link[sandwich:vcovCL]{sandwich::meatCL()}} if
\code{cluster} is defined, if \code{cluster = NULL}, they are passed into
\code{\link[sandwich:meat]{sandwich::meat()}}.}
}
\value{
An object of class \code{"chantrics"} inheriting from class \code{"chandwich"}.
See the documentation provided with \code{\link[chandwich:adjust_loglik]{chandwich::adjust_loglik()}}.
}
\description{
This function adjusts the loglikelihood of fitted model objects based on
Chandler and Bate (2007). It is a
generic function for different types of models, which are listed in
\strong{Supported models}. This section also contains links to function-specific
help pages.
}
\details{
If \code{use_vcov = TRUE}, the current default, the function will test
whether a \code{vcov} S3 method exists for \code{x}, and will take the
variance-covariance matrix from there. Otherwise, or if \code{use_vcov = FALSE}
the variance-covariance matrix of the MLE is estimated inside
\code{\link[chandwich:adjust_loglik]{chandwich::adjust_loglik()}} using \code{\link[stats:optim]{stats::optimHess()}}.
}
\section{Supported models}{

\itemize{
\item \link{glm}
\item \link{glm.nb}
}
}

\section{Available methods}{


\code{"chantrics"} objects have the following methods available to them:
\itemize{
\item \code{\link[stats:AIC]{AIC}}
\item \code{\link[=anova.chantrics]{anova}}
\item \code{\link[=alrtest]{alrtest}} - Adjusted Likelihood ratio tests
\item \code{\link[chandwich:coef.chandwich]{coef}}
\item \code{\link[lmtest:coeftest]{lmtest::coeftest}} - \(z\) tests for all
coefficients
\item \code{\link[chandwich:confint.chandwich]{confint}}
and \code{\link[chandwich:plot.confint]{plot.confint}} - confidence intervals for
all coefficients, and diagnostics plots for \code{confint()}.
\item \code{\link[chandwich:conf_intervals]{conf_intervals}} - enhanced confidence
interval reports
\item \code{\link[chandwich:conf_region]{conf_region}} - two-dimensional confidence
regions
\item \code{\link[stats:df.residual]{df.residual}}
\item \code{\link[stats:fitted.values]{fitted}}
\item \code{\link[chandwich:logLik.chandwich]{logLik}}
\item \code{\link[=logLik_vec]{logLik_vec}}
\item \code{\link[stats:nobs]{nobs}}
\item \code{\link[chandwich:plot.chandwich]{plot}}
\item \code{\link[=predict]{predict}}
\item \code{\link[chandwich:print.chandwich]{print}}
\item \code{\link[=residuals.chantrics]{residuals}}
\item \code{\link[chandwich:summary.chandwich]{summary}}
\item \code{\link[stats:terms]{terms}}
\item \code{\link[chandwich:vcov.chandwich]{vcov}}
}
}

\section{Examples}{
 See the model-specific pages in the \emph{supported models}
section.
}

\references{
R. E. Chandler and S. Bate, Inference for clustered data using the
independence loglikelihood, Biometrika, 94 (2007), pp. 167–183.
\doi{10.1093/biomet/asm015}.
}
\seealso{
\code{\link[lax:alogLik]{lax::alogLik()}} supports adjustment for user-supplied objects.
}
