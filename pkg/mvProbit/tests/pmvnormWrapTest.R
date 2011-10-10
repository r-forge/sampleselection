library( "mvProbit" )
library( "miscTools" )

# covariance matrix
sigma <- symMatrix( c( 1, 0.2, 0.4, 1, -0.1, 1 ) )

######## only upper ##########
upper <- c( -0.3, 0.7, -0.5 )
# Genz + Bretz (default)
pug <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = GenzBretz(), random.seed = 123 )
print( pug )

# Miwa (as function)
pum <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( pum )
all.equal( pug, pum, check.attributes = FALSE )

# Miwa (as object returned from function Miwa())
pum1 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = Miwa(), random.seed = 123 )
all.equal( pum, pum1 )

# Miwa (as character string)
pum2 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "Miwa", random.seed = 123 )
all.equal( pum, pum2 )

# TVPACK
put <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = TVPACK, random.seed = 123 )
print( put )
all.equal( pug, put, check.attributes = FALSE )
all.equal( pum, put, check.attributes = FALSE )

# GHK
pughk <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "ghk", random.seed = 123, nGHK = 1000 )
print( pughk )
all.equal( pug, pughk, check.attributes = FALSE )
all.equal( pum, pughk, check.attributes = FALSE )

# GHK, lower precision
pughk1 <- mvProbit:::pmvnormWrap( upper = upper, sigma = sigma, 
   algorithm = "ghk", random.seed = 123, nGHK = 100 )
print( pughk1 )
all.equal( pughk, pughk1, check.attributes = FALSE )
all.equal( pug, pughk1, check.attributes = FALSE )


######## only lower ##########
lower <- c( -0.7, 0.3, -0.9 )
# Genz + Bretz (default)
plg <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = GenzBretz, random.seed = 123 )
print( plg )

# Miwa
plm <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( plm )
all.equal( plg, plm, check.attributes = FALSE )

# TVPACK
plt <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = TVPACK, random.seed = 123 )
print( plt )
all.equal( plg, plt, check.attributes = FALSE )
all.equal( plm, plt, check.attributes = FALSE, tolerance = 1e-10 )

# GHK
plghk <- mvProbit:::pmvnormWrap( lower = lower, sigma = sigma, 
   algorithm = "GHK", random.seed = 123, nGHK = 1000 )
print( plghk )
all.equal( plg, plghk, check.attributes = FALSE )
all.equal( plm, plghk, check.attributes = FALSE )


######## partly lower, partly upper ##########
lower2 <- c( -Inf, 0.3, -Inf )
upper2 <- c( -0.3, Inf, -0.5 )
# Genz + Bretz (default)
pbg <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = GenzBretz(), random.seed = 123 )
print( pbg )

# Miwa
pbm <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = Miwa, random.seed = 123 )
print( pbm )
all.equal( pbg, pbm, check.attributes = FALSE )

# GHK
pbghk <- mvProbit:::pmvnormWrap( lower = lower2, upper = upper2, sigma = sigma, 
   algorithm = "GHK", random.seed = 123, nGHK = 1000 )
print( pbghk )
all.equal( pbg, pbghk, check.attributes = FALSE )
all.equal( pbm, pbghk, check.attributes = FALSE )


