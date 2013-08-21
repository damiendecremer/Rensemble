crpss <- function(data, crpssfile) {
  # Computes continuous ranked probability skill score (CRPSS) 
  # for a collection of ensemble forecasts
  # as proposed in Goddard et al (2013)
  #
  # ensembles are transformed to Gaussian forecast distributions
  # using sample mean and sample stdev
  #  
  # Usage: crpss(data,crpssfile)
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
  #    crpssfile: string of charaters with the name of the file where
  #             the CRPSS will be written 
  #  
  # Author: Stefan Siegert <s.siegert@exeter.ac.uk>
  
  data[data < -999] <- NA

  x <- data[, 3] # extracts observation column
  y <- data[, 4:(dim(data)[2])] # extracts ensemble forecasts
  ensmean <- apply(y, 1, mean, na.rm=TRUE)
  enssd <- apply(y, 1, sd, na.rm=TRUE)

  gauss.crps <- function(mu, sigma, obs) {
    sigma * (1/sqrt(pi) - 2 * dnorm(obs, mean=mu, sd=sigma) - 
    (obs - mu) / sigma * ( 2 * pnorm(obs, mean=mu, sd=sigma) - 1))
  }

  # calculate in-sample anomalies 
  x.anom <- x - mean(x, na.rm=TRUE)
  ensmean.anom <- ensmean - mean(ensmean)

  # remove conditional bias
  s.x <- sd(x.anom, na.rm=TRUE)
  s.y <- sd(ensmean.anom, na.rm=TRUE)
  r.xy <- cor(x.anom, ensmean.anom, use="pairwise.complete.obs")
  ensmean.anom.deb <- r.xy * s.x / s.y * ensmean.anom

  # calculate ensemble crps
  crps.e <- mean(gauss.crps(ensmean.anom.deb, enssd, x), na.rm=TRUE)
  crps.c <- mean(gauss.crps(0, s.x, x), na.rm=TRUE)

  #calculate CRPSS
  if (crps.c == 0) {
      if (crps.e == 0) {
        crpss <- 0 # i.e. ensemble mean and clim. are equally good
      } else {
        crpss <- -Inf
      }
  } else {
    crpss <- 1 - crps.e / crps.c
  }

  # write output file and return
  out <- crpss
  write.table(out, file=crpssfile)
  list(out=out)
} 

