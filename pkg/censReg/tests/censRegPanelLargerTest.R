library( "censReg" )
library( "plm" )

options( digits = 5 )

nId <- 100
nTime <- 3

set.seed( 123 )
pData <- data.frame(
   id = rep( paste( "F", 1:nId, sep = "_" ), each = nTime ),
   time = rep( 1980 + 1:nTime, nId ) )
pData$ui <- rep( rnorm( nId ), each = nTime )
pData$x1 <- rnorm( nId * nTime )
pData$x2 <- runif( nId * nTime )
pData$ys <- -1 + pData$ui + 2 * pData$x1 + 3 * pData$x2 + rnorm( nId * nTime )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
pData <- pdata.frame( pData, c( "id", "time" ), stringsAsFactors = FALSE )


# ## Newton-Raphson method
# randEff <- censReg( y ~ x1 + x2, data = pData )
# maxLik:::summary.maxLik( randEff )

## BHHH method
randEffBhhh <- censReg( y ~ x1 + x2, data = pData, method = "BHHH" )
print( maxLik:::summary.maxLik( randEffBhhh ), digits = 2 )

## BFGS method (optim)
randEffBfgs <- censReg( y ~ x1 + x2, data = pData, method = "BFGS" )
print( maxLik:::summary.maxLik( randEffBfgs ), digits = 2 )

## BFGS method (R)
randEffBfgsr <- censReg( y ~ x1 + x2, data = pData, method = "BFGSR" )
print( maxLik:::summary.maxLik( randEffBfgsr ), digits = 2 )

