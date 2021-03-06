predict.snqProfitEst <- function( object, newdata = object$data,
   se.fit = FALSE, se.pred = FALSE, interval = "none", level = 0.95, ... ) {

   nNetput <- length( object$pMeans )
   nFixed  <- length( object$fMeans )
   nObsOld <- nrow( object$data )
   nObsNew <- nrow( newdata )

   modelData <- .snqProfitModelData( data = newdata,
      weights = object$weights, pNames = object$pNames,
      qNames = object$qNames, fNames = object$fNames,
      ivNames = object$ivNames, form = object$form,
      netputScale = object$scalingFactors )
   system <- snqProfitSystem( nNetput, nFixed, form = object$form,
      profit = TRUE )
   restrict <- snqProfitRestrict( nNetput, nFixed, object$form )
   nCoefPerEq <- nrow( restrict ) / nNetput

   x <- list()
   result <- data.frame( obsNo = 1:nObsNew )
   for( i in 1:nNetput ) {
      x[[ i ]] <- model.matrix( system[[ i ]], modelData ) %*%
         restrict[ ( ( i - 1 ) * nCoefPerEq + 1 ):( i * nCoefPerEq ), ]
      result[[ object$qNames[ i ] ]] <- c( x[[ i ]] %*%
          object$coef$liCoef )
      if( se.fit || interval == "confidence" ) {
         result[[ paste( object$qNames[ i ], ".se.fit", sep = "" ) ]] <-
            diag( x[[ i ]] %*% object$coef$liCoefCov %*%
               t( x[[ i ]] ) )^0.5
      }
      if( se.pred || interval == "prediction" ) {
         result[[ paste( object$qNames[ i ], ".se.pred", sep = "" ) ]] <-
            diag( x[[ i ]] %*% object$coef$liCoefCov %*%
               t( x[[ i ]] ) + object$est$rcov[ i, i ] )^0.5
      }
      if( interval != "none" ) {
         if( object$est$probdfsys ) {
            tval   <- qt( 1 - ( 1- level )/2, object$est$df )
         } else {
            tval   <- qt( 1 - ( 1- level )/2, object$est$eq[[i]]$df )
         }
         if( interval == "confidence" ) {
            seName <- paste( object$qNames[ i ], ".se.fit", sep = "" )
         } else if( interval == "prediction" ) {
            seName <- paste( object$qNames[ i ], ".se.pred", sep = "" )
         } else {
            stop( "argument 'interval' must be either 'none', 'confidence'",
               " or 'prediction'" )
         }
         result[[ paste( object$qNames[ i ], ".lwr", sep="" ) ]] <-
            result[[ object$qNames[ i ] ]] - ( tval * result[[ seName ]] )
         result[[ paste( object$qNames[ i ], ".upr", sep="" ) ]] <-
            result[[ object$qNames[ i ] ]] + ( tval * result[[ seName ]] )
         if( !se.fit && interval == "confidence" ) result[[ seName ]] <- NULL
         if( !se.pred && interval == "prediction" ) result[[ seName ]] <- NULL
      }
   }
   if( !( "obsNo" %in% object$qNames ) ) result$obsNo <- NULL

   i <- nNetput + 1
   x[[ i ]] <- model.matrix( system[[ i ]], modelData )
   result$profit <- c( x[[ i ]] %*% object$coef$allCoef )
   if( se.fit || interval == "confidence" ) {
      result[[ "profit.se.fit" ]] <- diag( x[[ i ]] %*%
         object$coef$allCoefCov %*% t( x[[ i ]] ) )^0.5
   }
   if( se.pred || interval == "prediction" ) {
      s2 <- sum( residuals( object )$profit^2 ) / nObsOld
      result[[ "profit.se.pred" ]] <- diag( x[[ i ]] %*%
         object$coef$allCoefCov %*% t( x[[ i ]] ) +
         s2 )^0.5
   }
   if( interval != "none" ) {
      if( object$est$probdfsys ) {
         tval   <- qt( 1 - ( 1- level )/2, object$est$df )
      } else {
         tval   <- qt( 1 - ( 1- level )/2, nObsOld )
      }
      if( interval == "confidence" ) {
         seName <- "profit.se.fit"
      } else if( interval == "prediction" ) {
         seName <- "profit.se.pred"
      } else {
         stop( "argument 'interval' must be either 'none', 'confidence'",
            " or 'prediction'" )
      }
      result[[ "profit.lwr" ]] <- result[[ "profit" ]] -
         ( tval * result[[ seName ]] )
      result[[ "profit.upr" ]] <- result[[ "profit" ]] +
         ( tval * result[[ seName ]] )
      if( !se.fit && interval == "confidence" ) result[[ seName ]] <- NULL
      if( !se.pred && interval == "prediction" ) result[[ seName ]] <- NULL
   }

   return( result )
}

predict.snqProfitImposeConvexity <- function( object, newdata = object$data,
   se.fit = FALSE, se.pred = FALSE, interval = "none", level = 0.95, ... ) {

   if( is.null( object$sim ) ) {
      if( se.fit ) {
         warning( "setting argument 'se.fit' to 'FALSE' because",
            " standard errors are not available" )
         se.fit <- FALSE
      }
      if( se.pred ) {
         warning( "setting argument 'se.pred' to 'FALSE' because",
            " standard errors are not available" )
         se.pred <- FALSE
      }
      if( interval != "none" ) {
         warning( "setting argument 'interval' to 'none' because",
            " standard errors are not available" )
         interval <- "none"
      }
   }

   result <- predict.snqProfitEst( object, newdata = newdata,
      se.fit = se.fit, se.pred = se.pred, interval = interval,
      level = level, ... )

   return( result )
}
