reliability_new1 <- function(data, u, nbins=10, nboot=1000, 
                             threshold = TRUE, reliabfile,
                             graphvaluefile, maintitle="") 
{
  #
  # Plot reliability diagram for an ensemble forecast
  #
  # Usage: reliability_new1(data, u, nbins, nboot, threshold, reliabfile, graphvaluefile, maintitle)
  #
  # Arguments:
  #    data: a table with the following structure
  #          first  column = year
  #	     second column = month
  #	     third  column = observation
  #	     forth  column = ensemble member 1
  #	     fith   column = ensemble member 2
  #	     .
  #	     .
  #	     .
  #	     last   column = ensemble member k
  #
  #    u: threshold value or quantile (e.g. u=50 for the 0.5 quantile, 
  #       i.e. the median of the observations) that defined the event of 
  #       interest to be forecast
  #
  #    nbins: number of probability bins 
  #
  #    nboot: number of bootstrap cycles to estimate consistency bars
  #
  #    threshold: Logical. If TRUE (Default) uses the value u as threshold, if
  #               FALSE u should have the value of the quantile to be computed
  #  
  #
  #    reliabfile: string of charaters with the name of the file where
  #            the threshold value or quantile u will be written
  #
  # graphvaluefile: string of charaters with the name of the file where
  #                 the data used to produce reliability diagram plot will
  #                 be written
  #                   
  #    maintitle: String containing the text for the reliability diagram title
  #
  #
  # change log:
  #
  #  2013/08/20:
  #  * points are plotted at in-bin-averages, not at bin centres
  #  * legend has been removed
  #  * consistency bars have been added, calculated by a resampling technique
  #  * see Broecker (2007) http://dx.doi.org/10.1175/WAF993.1 for details
  #  * the bars are pointwise 2.5% ... 97.5% intervals around the hypothesis of reliability
  #  * dependency on package "verification" was removed
  #
  # Author: Stefan Siegert <s.siegert@exeter.ac.uk>
  #
  # based on previous version by Caio Coelho and the routine 
  # verification::reliability.plot.default
  #


  data[data == -999.9] <- NA

  x <- data[, 3] # observations
  y <- data[, 4:(dim(data)[2])] # ensemble forecasts

  # delete ensemble columns which have only NA
  na.inds <- which(apply(y, 2, function(z) all(is.na(z))))
  if (length(na.inds) > 0) {
    y <- y[, -na.inds]
  }
  # delete rows which contain at least one NA
  na.inds <- which(is.na(rowSums(cbind(x, y))))
  if (length(na.inds) > 0) {
    x <- x[-na.inds]
    y <- y[-na.inds, ]
  }

  # threshold or quantile?
  if(!threshold) u <- quantile(x, u / 100)

  # probability forecasts and binary observations
  probfcsts <- rowMeans(y <= u)
  binobs <- 1 * (x <= u)
  n <- length(binobs)

  # estimate refinement function
  brx <- seq(0, 1, length.out=nbins+1) + 
         c(-.Machine$double.eps, rep(0, nbins-1), .Machine$double.eps)
  h <- hist(probfcsts, breaks=brx, plot=FALSE)$counts        

  # estimate calibration function
  g <- hist(probfcsts[binobs==1], breaks=brx, plot=FALSE)$counts
  obar.i <- g / h 
  
  # calculate in-bin averages
  p.bins <- as.numeric(cut(probfcsts, breaks=brx))
  p.avgs <- sapply(seq(nbins), 
                   function(ii) mean(probfcsts[p.bins == ii], na.rm=TRUE))

  # consistency resampling (broecker and smith 2007)
  resamp.mat <- matrix(nrow=0, ncol=nbins)
  for (i in 1:nboot) {
    p.hat <- sample(x=probfcsts, size=n, replace=TRUE)
    x.hat <- rbinom(n=n, size=1, prob=p.hat)
    hh <- hist(p.hat, breaks=brx, plot=FALSE)$counts        
    gg <- hist(p.hat[x.hat==1], breaks=brx, plot=FALSE)$counts
    resamp.mat <- rbind(resamp.mat, gg / hh)
  }
  cons.bars <- apply(resamp.mat, 2, 
                     function(z) quantile(z, c(.025, .975), na.rm=TRUE))

  # reliability plot
  old.par <- par(no.readonly = TRUE) 
  on.exit(par(old.par))
  plot(NULL, xlim = c(0,1), ylim = c(0,1),
     xlab= expression(paste("Forecast probability, ", y[i])),
     ylab=expression(paste("Observed relative frequency, ", bar(o)[1])),
     main=maintitle)
  # consistency bars
  for (i in 1:length(p.avgs)) {
      lines(c(p.avgs[i], p.avgs[i]), cons.bars[, i], col="#CCCCCC", lwd=6)
  }

  # reliability points and diagonal
  points(p.avgs, obar.i, type = "b", col = "red", lty = 1, lwd = 2)
  abline(0,1)
  
  # refinement histogram in lower corner
  pp<- par("plt")
  par("plt" = c(pp[2] - 0.2 , pp[2],  pp[3], pp[3]+ 0.2) )
  par(new = TRUE)
  barplot(h, axes = FALSE, axisnames = FALSE)
  #hist(probfcsts, axes=FALSE, main=NA, xlab=NA, ylab=NA, breaks=brx, col="gray")
  #lines(density(probfcsts, width=0.01), col="black")
  axis(4)
  box() 
  
  # save files
  out<-cbind(p.avgs,obar.i,h,t(cons.bars))
  write.table(round(u,2), file=reliabfile)
  write.table(round(out,2),file=graphvaluefile)
}
