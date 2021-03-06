\name{censReg}
\alias{censReg}
\alias{print.censReg}
\title{Censored Regression (Tobit) Model}
\description{
   Fitting a model with a censored dependent variable.
}

\usage{
censReg( formula, left = 0, right = Inf, data = sys.frame( sys.parent()),
   subset = NULL, start = NULL, nGHQ = 8, logLikOnly = FALSE, ... )

\method{print}{censReg}( x, logSigma = TRUE, digits = 4, ... )
}

\arguments{
   \item{formula}{an object of class \code{"formula"}:
      a symbolic description of the model to be fitted.}
   \item{left}{left limit for the censored dependent variable;
      if set to \code{-Inf}, the dependent variable is assumed to be
      not left-censored; defaults to zero (classical Tobit model).}
   \item{right}{right limit for the censored dependent variable;
      if set to \code{Inf}, the dependent variable is assumed to be
      not right-censored; defaults to \code{Inf} (classical Tobit model).}
   \item{data}{an optional data frame.
      If argument \code{data} is of class \code{"pdata.frame"},
      a panel-model is estimated.}
   \item{subset}{an optional vector specifying a subset of observations 
      to be used in the fitting process.}
   \item{start}{an optional vector of initial parameters for the ML estimation
      (intercept, slope parameters, logarithm of the standard deviation
      of the individual effects (only for random-effects panel data models), 
      and logarithm of the standard deviation of the general disturbance term);
      if \code{start} is not specified, initial values are taken from an OLS
      estimation or (uncensored) random-effects panel data estimation.}
   \item{nGHQ}{number of points used in the Gauss-Hermite quadrature,
      which is used to compute the log-likelihood value in case of the
      random effects model for panel data.}
   \item{logLikOnly}{logical. If \code{TRUE}, the model is not estimated
      but the log-likelihood contributions of all observations/individuals
      and the corresponding gradients 
      calculated at the parameters specified by argument \code{start}
      are returned.}
   \item{x}{object of class \code{censReg} (returned by \code{censReg}).}
   \item{logSigma}{logical value indicating whether the variance(s)
      of the model should be printed logarithmized
      (see \code{\link{coef.censReg}}).}
   \item{digits}{positive integer specifiying the minimum number of
      significant digits to be printed
      (see \code{\link[base]{print.default}}).}
   \item{\dots}{additional arguments for \code{censReg} are passed to
      \code{\link[maxLik]{maxLik}};
      additional arguments for \code{print.censReg} are currently ignored.}
}

\value{
   If argument \code{logLikOnly} is \code{FALSE} (default),
   \code{censReg} returns an object of class \code{"censReg"}
   inheriting from class \code{"maxLik"}.
   The returned object contains the same components as objects
   returned by \code{\link[maxLik]{maxLik}} and additionally
   the following components:
   \item{call}{the matched call.}
   \item{terms}{the model terms.}
   \item{nObs}{a vector containing 4 integer values:
      the total number of observations,
      the number of left-censored observations,
      the number of uncensored observations, and
      the number of right-censored observations.}
   \item{df.residual}{degrees of freedom of the residuals.}
   \item{start}{vector of starting values.}
   \item{left}{left limit of the censored dependent variable.}
   \item{right}{right limit of the censored dependent variable.}
   \item{xMean}{vector of mean values of the explanatory variables.}

   In contrast, 
   if argument \code{logLikOnly} is \code{TRUE},
   \code{censReg} returns a vector 
   of the log-likelihood contributions of all observations/individuals.
   This vector has an attribute \code{"gradient"},
   which is a matrix containing the gradients of the log-likelihood contributions
   with respect to the parameters.
}

\details{
   The model is estimated by Maximum Likelihood (ML)
   assuming a Gaussian (normal) distribution of the error term.
   The maximization of the likelihood function is done
   by function \code{\link[maxLik]{maxLik}} of the \pkg{maxLik} package.
   An additional argument \code{method} can be used to specify
   the optimization method used by \code{\link[maxLik]{maxLik}},
   e.g.\ \code{"Newton-Raphson"}, \code{"BHHH"}, \code{"BFGS"},
   \code{"SANN"} (for simulated annealing), or
   \code{"NM"} (for Nelder-Mead).
}

\references{
   Greene, W.H. (2008):
   \emph{Econometric Analysis}, Sixth Edition, Prentice Hall, p. 871-875.

   Kleiber, C. and Zeileis, A. (2008):
   \emph{Applied Econometrics with R}, Springer, p. 141-143.

   Tobin, J. (1958):
   Estimation of Relationships for Limited Dependent Variables.
   \emph{Econometrica} 26, p. 24-36.
}

\note{
   When the censored regression model is estimated,
   the log-likelihood function is maximized with respect
   to the coefficients and the \emph{logarithm(s)}
   of the variance(s).
}

\author{Arne Henningsen}

\seealso{\code{\link{summary.censReg}}, \code{\link{coef.censReg}},
   \code{\link[AER]{tobit}}, \code{\link[sampleSelection]{selection}}}

\examples{
## Kleiber & Zeileis ( 2008 ), page 142
data( "Affairs", package = "AER" )
estResult <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs )
print( estResult )

## Kleiber & Zeileis ( 2008 ), page 143
estResultBoth <- censReg( affairs ~ age + yearsmarried + religiousness +
   occupation + rating, data = Affairs, right = 4 )
print( estResultBoth )
}

\keyword{models}
\keyword{regression}
