mvProbitExp <- function( formula, coef, sigma, data, yNames = NULL, 
   cond = FALSE, random.seed = 123, ... ) {

   # checking argument 'cond'
   if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be logical" )
   } else if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical values" )
   }

   # checking argument 'random.seed'
   if( !is.numeric( random.seed ) ) {
      stop( "argument 'random.seed' must be numerical" )
   } else if( length( random.seed ) != 1 ) {
      stop( "argument 'random.seed' must be a single numerical values" )
   }

   # checking argument 'formula'
   if( is.list( formula ) ) {
      stop( "using different regressors for the dependent variables",
         " has not been implemented yet. Sorry!" )
   } else if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a data frame" )
   }

   # preparing model matrix
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )

   # checking argument 'sigma'
   if( !is.matrix( sigma ) ) {
      stop( "argument 'sigma' must be a matrix" )
   } else if( nrow( sigma ) != ncol( sigma ) ) {
      stop( "argument 'sigma' must be a quadratic matrix" )
   } else if( !isSymmetric( sigma ) ) {
      stop( "argument 'sigma' must be a symmetric matrix" )
   }

   # number of dependent variables
   nDep <- ncol( sigma )

   # number of regressors
   nReg <- ncol( xMat )

   # number of coefficients
   nCoef <- nDep * nReg

   # number of observations
   nObs <- nrow( xMat )

   # checking argument 'coef'
   if( !is.vector( coef, mode = "numeric" ) ) {
      stop( "argument 'coef' must be a numeric vector" )
   } else if( length( coef ) != nCoef ) {
      stop( "argument coef must have ", nCoef, " elements" )
   }

   # separating coefficients for different equations
   betaEq <- list()
   for( i in 1:nDep ) {
      betaEq[[ i ]] <- coef[ ( ( i - 1 ) * nReg + 1 ):( i * nReg ) ]
   }

   # checking argument 'yNames'
   if( !is.null( yNames ) ) {
      if( length( yNames ) != nDep ) {
         stop( "argument 'yNames' must be either 'NULL'",
            " or must have the same length as the number of rows and columns",
            " of argument 'sigma'" )
      }
   }

   # preparing matrix of responses
   if( !is.null( yNames ) ) {
      yMat <- matrix( NA, nrow = nObs, ncol = nDep )
      for( k in 1:nDep ) {
         yMat[ , k ] = data[[ yNames[ k ] ]]
      }
   }

   # calculating linear predictors
   xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
   for( i in 1:nDep ) {
      xBeta[ , i ] <- xMat %*% betaEq[[ i ]]
   }

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by pmvnorm)
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   if( cond ) {
      # conditional expectations
      result <- matrix( NA, nrow = nObs, ncol = nDep )
      if( is.null( yNames ) ) {
         # assuming that all other dependent variables are one
         for( i in 1:nObs ) {
            for( k in 1:nDep ) {
               result[ i, k ] <- 
                  pmvnorm( upper = xBeta[ i, ], sigma = sigma, ... ) / 
                  pmvnorm( upper = xBeta[ i, -k ], sigma = sigma[ -k, -k ], ... )
            }
         }
      } else {
         # assuming that all other dependent variables are as observed
         for( i in 1:nObs ){
            for( k in 1:nDep ) {
               ySign <- 2 * yMat[ i, ] - 1
               ySign[ k ] <- 1
               xBetaTmp <- xBeta[ i, ] * ySign
               sigmaTmp <- diag( ySign ) %*% sigma %*% diag( ySign )
               result[ i, k ] <- 
                  pmvnorm( upper = xBetaTmp, sigma = sigmaTmp, ... ) / 
                  pmvnorm( upper = xBetaTmp[ -k ], sigma = sigmaTmp[ -k, -k ], ... )
            }
         }
      }
   } else {
      result <- pnorm( xBeta )
   }

   if( !is.null( yNames ) ) {
      colnames( result ) <- yNames
   }

   result <- as.data.frame( result )

   return( result )
}