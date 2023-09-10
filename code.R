################################################################################
###
### Code illustrating the analyses in:
###
### The effect of temporal data aggregation to assess the impact of changing temperatures in Europe: an epidemiological modelling study
### The Lancet Regional Health - Europe
###
### Contact: Joan Ballester (joan.ballester@isglobal.org)
###
### Note: This is a simplified code with sample data, results are expected to
###       differ from those published in the article
### 
################################################################################

rm( list = ls() );
cat("\014");

# Required Libraries
suppressMessages( library(dlnm) ); # crossbasis
suppressMessages( library(splines) ); # ns, bs



################################################################################
### Parameter Definition
################################################################################

# Exposure-Response: Temperature Knot Percentiles
VAR_PRC = c(10,75,90) / 100;

# Lag-Response: Minimum and Maximum Lags (in Weeks)
MIN_LAG =  0; if( 0 > MIN_LAG ){ stop( "ERROR: Invalid MIN_LAG !!!" ); }
MAX_LAG = 28; if( MIN_LAG >= MAX_LAG | MAX_LAG %% 7 != 0 | MAX_LAG %% 4 != 0 ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }

# Degrees of Freedom per Year for the Seasonal and Long-Term Trends
DF_SEAS = 8; if( DF_SEAS <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }



################################################################################
### Data Preparation
################################################################################

# Reading the Original Chicago Data
DATA_ORI = chicagoNMMAPS;
DATA_ORI$mort = DATA_ORI$death;
DATA_ORI$death = NULL;

# Identifying the Range of Dates with Full Monday-To-Sunday Weeks
iTIM1 =      which( DATA_ORI$dow == "Monday" )  [1];
iTIM2 = rev( which( DATA_ORI$dow == "Sunday" ) )[1];
if( ( iTIM2 - iTIM1 + 1 ) %% 7 != 0 ){ stop("ERROR: Invalid Time Vector I/II !!!"); }

# Identifying the Range of Dates with Full Monday-To-Sunday MAX_LAG-Periods
iTIM2 = iTIM2 - ( iTIM2 - iTIM1 + 1 ) %% MAX_LAG;
if( ( iTIM2 - iTIM1 + 1 ) %% MAX_LAG != 0 ){ stop("ERROR: Invalid Time Vector II/II !!!"); }

# Data Frames: Daily, Weekly, 2-Weekly and 4-Weekly
DATA_AGG = vector( "list", 4 );
names( DATA_AGG ) = c( "Daily", "Weekly", "2-Weekly", "4-Weekly" );

# Creating the Daily Data Table
DATA_AGG[[1]] = DATA_ORI[ iTIM1:iTIM2, c( "date", "dow", "temp", "mort" ) ];
rm(DATA_ORI);

# Creating the Weekly Data Table
nTIM_AGG2 = ( iTIM2 - iTIM1 + 1 ) / 7;
vTIM_AGG2 = rep( seq( nTIM_AGG2 ), each = 7 );
DATA_AGG[[2]] = data.frame( time =                  1:nTIM_AGG2                 ,   # Time ID
                            temp = tapply( DATA_AGG[[1]]$temp, vTIM_AGG2, mean ),   # Temperature
                            mort = tapply( DATA_AGG[[1]]$mort, vTIM_AGG2, sum  ) ); # Mortality
rm(nTIM_AGG2,vTIM_AGG2);

# Creating the 2-Weekly Data Table
nTIM_AGG3 = ( iTIM2 - iTIM1 + 1 ) / 14;
vTIM_AGG3 = rep( seq( nTIM_AGG3 ), each = 14 );
DATA_AGG[[3]] = data.frame( time =                  1:nTIM_AGG3                 ,   # Time ID
                            temp = tapply( DATA_AGG[[1]]$temp, vTIM_AGG3, mean ),   # Temperature
                            mort = tapply( DATA_AGG[[1]]$mort, vTIM_AGG3, sum  ) ); # Mortality
rm(nTIM_AGG3,vTIM_AGG3);

# Creating the 4-Weekly Data Table
nTIM_AGG4 = ( iTIM2 - iTIM1 + 1 ) / 28;
vTIM_AGG4 = rep( seq( nTIM_AGG4 ), each = 28 );
DATA_AGG[[4]] = data.frame( time =                  1:nTIM_AGG4                 ,   # Time ID
                            temp = tapply( DATA_AGG[[1]]$temp, vTIM_AGG4, mean ),   # Temperature
                            mort = tapply( DATA_AGG[[1]]$mort, vTIM_AGG4, sum  ) ); # Mortality
rm(nTIM_AGG4,vTIM_AGG4);

rm(iTIM1,iTIM2);



################################################################################
### Location-Specific Associations
################################################################################

# Model Cross-Predictions
CROSS_PRED = vector( "list", length(DATA_AGG) );
names(CROSS_PRED) = names(DATA_AGG);

for( iAGG in 1:length(DATA_AGG) ){
  
  # Model Formula and Temperature Cross-Basis
  if( names(DATA_AGG)[iAGG] == "Daily" ){
    
    FORMULA = mort ~ ns( date, df = round( DF_SEAS * length(date) * 1 / 365.25 ) ) + CROSS_BASIS + dow;
    
    CROSS_BASIS = crossbasis( DATA_AGG[[iAGG]]$temp,
                              c( MIN_LAG, MAX_LAG ),
                              argvar = list( fun = "ns", knots = quantile( DATA_AGG[[iAGG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( knots = logknots( c( MIN_LAG, MAX_LAG ), 3 ) ) );
    
  }else if( names(DATA_AGG)[iAGG] == "Weekly" ){
    
    FORMULA = mort ~ ns( time, df = round( DF_SEAS * length(time) * 7 / 365.25 ) ) + CROSS_BASIS;
    
    CROSS_BASIS = crossbasis( DATA_AGG[[iAGG]]$temp,
                              c( MIN_LAG, MAX_LAG / 7 ),
                              argvar = list( fun = "ns", knots = quantile( DATA_AGG[[iAGG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
    
  }else if( names(DATA_AGG)[iAGG] == "2-Weekly" ){
    
    FORMULA = mort ~ ns( time, df = round( DF_SEAS * length(time) * 14 / 365.25 ) ) + CROSS_BASIS;
    
    CROSS_BASIS = crossbasis( DATA_AGG[[iAGG]]$temp,
                              c( MIN_LAG, MAX_LAG / 14 ),
                              argvar = list( fun = "ns", knots = quantile( DATA_AGG[[iAGG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
    
  }else if( names(DATA_AGG)[iAGG] == "4-Weekly" ){
    
    FORMULA = mort ~ ns( time, df = round( DF_SEAS * length(time) * 28 / 365.25 ) ) + CROSS_BASIS;
    
    CROSS_BASIS = crossbasis( DATA_AGG[[iAGG]]$temp,
                              c( MIN_LAG, MAX_LAG / 28 ),
                              argvar = list( fun = "ns", knots = quantile( DATA_AGG[[iAGG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
    
  }else{
    stop( "ERROR: Invalid iAGG !!!" );
  }
  
  # Fitting the Model
  GLM_MODEL = glm( formula = FORMULA, DATA_AGG[[iAGG]], family = quasipoisson, na.action = "na.exclude" );
  rm(FORMULA);
  
  # Cumulative Exposure-Response without Centering
  suppressMessages( CROSS_PRED[[iAGG]] <- crosspred( CROSS_BASIS, GLM_MODEL, at = quantile( DATA_AGG[[iAGG]]$temp, seq(0,1,0.001), na.rm = TRUE ) ) );
  
  # Minimum Mortality Temperature
  MMT = CROSS_PRED[[iAGG]]$predvar[ which.min( CROSS_PRED[[iAGG]]$allRRfit ) ];
  
  # Cumulative Exposure-Response with Centering at the Minimum Mortality Temperature
  CROSS_PRED[[iAGG]] = crosspred( CROSS_BASIS, GLM_MODEL, at = quantile( DATA_AGG[[iAGG]]$temp, seq(0,1,0.001), na.rm = TRUE ), cen = MMT );
  rm(CROSS_BASIS,GLM_MODEL, MMT);
  
}
rm(iAGG);



################################################################################
### Plots
################################################################################

# Cumulative Exposure-Response

pdf( "./plot_RR.pdf", width = 4, height = 4 );
par( mex = 0.8, mgp = c(2.5,1,0), las = 0 );

plot   (    CROSS_PRED[[1]]$predvar                                  ,    CROSS_PRED[[1]]$allRRfit                                    , col = rgb(0.00,0.00,0.00    ), lwd = 4, lty = 1, type = "l", main = "Chicago", xlab = expression( paste( "Temperature (", degree, "C)" ) ), ylab = "Relative Risk", xlim = range( c( CROSS_PRED[[1]]$predvar, CROSS_PRED[[2]]$predvar, CROSS_PRED[[3]]$predvar, CROSS_PRED[[4]]$predvar ) ), ylim = range( c( CROSS_PRED[[1]]$allRRfit, CROSS_PRED[[2]]$allRRfit, CROSS_PRED[[3]]$allRRfit, CROSS_PRED[[4]]$allRRfit ) ), axes = T );
polygon( c( CROSS_PRED[[1]]$predvar, rev( CROSS_PRED[[1]]$predvar ) ), c( CROSS_PRED[[1]]$allRRlow, rev( CROSS_PRED[[1]]$allRRhigh ) ), col = rgb(0.00,0.00,0.00,1/3), border = FALSE );
polygon( c( CROSS_PRED[[2]]$predvar, rev( CROSS_PRED[[2]]$predvar ) ), c( CROSS_PRED[[2]]$allRRlow, rev( CROSS_PRED[[2]]$allRRhigh ) ), col = rgb(1.00,0.00,0.00,1/3), border = FALSE );
polygon( c( CROSS_PRED[[3]]$predvar, rev( CROSS_PRED[[3]]$predvar ) ), c( CROSS_PRED[[3]]$allRRlow, rev( CROSS_PRED[[3]]$allRRhigh ) ), col = rgb(0.00,0.00,1.00,1/3), border = FALSE );
polygon( c( CROSS_PRED[[4]]$predvar, rev( CROSS_PRED[[4]]$predvar ) ), c( CROSS_PRED[[4]]$allRRlow, rev( CROSS_PRED[[4]]$allRRhigh ) ), col = rgb(0.50,0.50,0.50,1/3), border = FALSE );
lines  (    CROSS_PRED[[1]]$predvar                                  ,    CROSS_PRED[[1]]$allRRfit                                    , col = rgb(0.00,0.00,0.00    ), lwd = 4, lty = 1 );
lines  (    CROSS_PRED[[2]]$predvar                                  ,    CROSS_PRED[[2]]$allRRfit                                    , col = rgb(1.00,0.00,0.00    ), lwd = 4, lty = 1 );
lines  (    CROSS_PRED[[3]]$predvar                                  ,    CROSS_PRED[[3]]$allRRfit                                    , col = rgb(0.00,0.00,1.00    ), lwd = 4, lty = 1 );
lines  (    CROSS_PRED[[4]]$predvar                                  ,    CROSS_PRED[[4]]$allRRfit                                    , col = rgb(0.50,0.50,0.50    ), lwd = 4, lty = 1 );
abline ( h = 1, lwd = 1, lty = 1, col = rgb(0.00,0.00,0.00) );
legend( "topright", names(DATA_AGG), col = rgb( c(0.00,1.00,0.00,0.50),c(0.00,0.00,0.00,0.50),c(0.00,0.00,1.00,0.50) ), lwd = 4, lty = 1, cex = 1.000, box.lty = 0, bg = "transparent" );

dev.off();

# "We are a way for the cosmos to know itself" (Carl Sagan)
# "We've tried nothing, and we're all out of ideas" (The Simpsons)
