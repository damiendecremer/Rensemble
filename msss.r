msss <- function(data, msssfile) {

  # Computes mean squared skill score (MSSS) of the ensemble mean
  # as proposed in Goddard et al (2013)
  #  
  # Usage: msss(data,msssfile)
  #
  # Arguments:
  #    data: a table with the following structure
  #           first  column = year
  #	     second column = month
  #	     third  column = observation
  #	     forth  column = ensemble member 1
  #	     fith   column = ensemble member 2
  #	     .
  #	     .
  #	     .
  #	     last   column = ensemble member k
  #
  #    msssfile: string of charaters with the name of the file where
  #             the MSSS will be written 
  #  
  # Author: Stefan Siegert <s.siegert@exeter.ac.uk>
  
  data[data < -999] <- NA

  x <- data[, 3] # extracts observation column
  y <- data[, 4:(dim(data)[2])] # extracts ensemble forecasts
  ensmean <- apply(y, 1, mean, na.rm=TRUE)

  # calculate in-sample anomalies 
  ensmean.anom <- ensmean - mean(ensmean, na.rm=TRUE)
  x.anom <- x - mean(x, na.rm=TRUE)

  # calculate MSE of the ensemble mean 
  mse.e <- mean((ensmean.anom - x.anom)^2, na.rm=TRUE)
  # calculate MSE of climatology (clim. = 0 by construction)
  mse.c <- mean(x.anom^2, na.rm=TRUE) 

  #calculate MSSS
  if (mse.c == 0) {
      if (mse.e == 0) {
        msss <- 0 # i.e. ensemble mean and clim. are equally good
      } else {
        msss <- -Inf
      }
  } else {
    msss <- 1 - mse.e / mse.c
  }

  # write output file and return
  out <- msss
  write.table(out, file=msssfile)
  list(out=out)
} 

