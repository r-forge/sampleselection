
options( digits = 3 )
library( "sampleSelection" )
set.seed(0)
                           # Note: library() may change RNG state!

## Leeman Lucas (and many others): binary outcome

N <- 500
rho <- 0.7
suppressPackageStartupMessages( library( "mvtnorm" ) )
eps <- rmvnorm(N, c(0,0), matrix(c(1,rho,rho,1), 2, 2) )
simDat <- data.frame( xs = runif(N) )
simDat$ysX <- 3 * simDat$xs + eps[,1]
simDat$ys <- simDat$ysX > 0
simDat$xo <- runif(N)
simDat$yoX <- -1 + 2 * simDat$xo + eps[,2]
simDat$yo <- factor( (simDat$yoX > 0) * (simDat$ys > 0))
                           # binary outcome, only observable if ys>0
print(table(simDat$ys, simDat$yo, exclude=NULL))

# estimation with BHHH method
ss <- selection( ys ~ xs, yo ~ xo, data = simDat, steptol = 1e-12 )
print( ss )
summary( ss )
coef( ss )
coef( ss, part = "outcome" )
coef( summary( ss ) )
coef( summary( ss ), part = "outcome" )
stdEr( ss )
vcov( ss )
vcov( ss, part = "outcome" )
nobs( ss )
nObs( ss )
round( fitted( ss ), 3 )
all.equal( fitted( ss ), fitted( ss, part = "outcome" ) )
round( fitted( ss, part = "selection" ), 3 )
round( residuals( ss ), 3 )
all.equal( residuals( ss ), residuals( ss, part = "outcome" ) )
all.equal( residuals( ss ),
   residuals( ss, part = "outcome", type = "deviance"  ) )
round( residuals( ss, type = "pearson" ), 3 )
all.equal( residuals( ss, type = "pearson" ),
   residuals( ss, part = "outcome", type = "pearson" ) )
round( residuals( ss, type = "deviance" ), 3 )
all.equal( residuals( ss, type = "deviance" ),
   residuals( ss, part = "outcome", type = "deviance" ) )
all.equal( residuals( ss, part = "outcome", type = "response" ),
   ( simDat$yo == 1 ) - fitted( ss, part = "outcome" ) )
round( residuals( ss, part = "selection" ), 3 )
all.equal( residuals( ss, part = "selection" ),
   residuals( ss, part = "selection", type = "deviance" ) )
round( residuals( ss, part = "selection", type = "pearson" ), digits = 3 )
round( residuals( ss, part = "selection", type = "response" ), digits = 3 )
all.equal( residuals( ss, part = "selection", type = "response" ),
   simDat$ys - fitted( ss, part = "selection" ) )
model.matrix( ss )
all.equal( model.matrix( ss ), model.matrix( ss, part = "outcome" ) )
model.matrix( ss, part = "selection" )
model.frame( ss )
logLik( ss )

# estimation with BFGS method
ssBFGS <- selection( ys ~ xs, yo ~ xo, data = simDat, maxMethod = "BFGS" )
print( ssBFGS )
summary( ssBFGS )
all.equal( coef( ssBFGS ), coef( ss ), tol = 1e-2 )
all.equal( stdEr( ssBFGS ), stdEr( ss ), tol = 1e-1 )
all.equal( vcov( ssBFGS ), vcov( ss ), tol = 1e-1 )
nobs( ssBFGS )
nObs( ssBFGS )
all.equal( fitted( ss ), fitted( ssBFGS ), tol = 1e-2 )
all.equal( fitted( ss, part = "selection" ),
   fitted( ssBFGS, part = "selection" ), tol = 1e-2 )
all.equal( fitted( ssBFGS ), fitted( ssBFGS, part = "outcome" ) )
all.equal( residuals( ss ), residuals( ssBFGS ), tol = 1e-2 )
all.equal( residuals( ss, part = "selection" ),
   residuals( ssBFGS, part = "selection" ), tol = 1e-2 )
all.equal( residuals( ssBFGS ), residuals( ssBFGS, part = "outcome" ) )
all.equal( model.matrix( ss ), model.matrix( ssBFGS ) )
all.equal( model.matrix( ss, part = "selection" ),
   model.matrix( ssBFGS, part = "selection" ) )
all.equal( model.matrix( ssBFGS ), model.matrix( ssBFGS, part = "outcome" ) )
all.equal( logLik( ss ), logLik( ssBFGS ), tol = 1e-3 )
