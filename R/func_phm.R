########## functions for PCH #############################
library(tidyr)
library(dplyr)

# function for creating stratum-specific time intervals so that all intervals have the same number of events
create_intervals <- function(time, event, S, historical, n.intervals){
  
  df <- data.frame(time=time, event=event, S=S)
  df <- df[order(df[,"time"],decreasing=FALSE),]
  
  unique_strat <- sort(unique(S))
  
  # loop through the strata
  interval_list <- list()
  for(k in 1:length(unique_strat)){
    
    df_s <- df %>% filter(S==unique_strat[k])
    ev_time <- df_s$time[df_s$event==1]
    if(length(historical) != 0){
      for(i in 1:length(historical)){
        data <- historical[[i]]
        df_data <- data.frame(time=data[["time"]], event=data[["event"]], S=data[["S"]])
        df_h <- df_data %>% filter(S==unique_strat[k])
        ev_time <- c(ev_time, df_h$time[df_h$event==1])
      }
    }
    intervals <- quantile(ev_time, probs=seq(from=0,to=1,length.out=n.intervals[k]+1)) 
    intervals[1] <- 0
    intervals[n.intervals[k]+1] <- Inf
    
    interval_list[[k]] <- intervals
  }
  
  return(interval_list)
  
}



# for design: function for creating stratum-specific time intervals so that all intervals have the same number of events
create_intervals_design <- function(historical, n.intervals){
  
  df <- historical[[1]]
  S <- df[["S"]]
  
  unique_strat <- sort(unique(S))
  
  # loop through the strata
  interval_list <- list()
  for(k in 1:length(unique_strat)){
    
    ev_time <- NULL
    if(length(historical) != 0){
      for(i in 1:length(historical)){
        data <- historical[[i]]
        df_data <- data.frame(time=data[["time"]], event=data[["event"]], S=data[["S"]])
        df_h <- df_data %>% filter(S==unique_strat[k])
        ev_time <- c(ev_time, df_h$time[df_h$event==1])
      }
    }
    intervals <- quantile(ev_time, probs=seq(from=0,to=1,length.out=n.intervals[k]+1)) 
    intervals[1] <- 0
    intervals[n.intervals[k]+1] <- Inf
    
    interval_list[[k]] <- intervals
  }
  
  return(interval_list)
  
}


create_tables_new <- function(time, event, X, S, change.points){
  
  X <- as.matrix(X)
  P <- ncol(X)
  dffull <- data.frame(time=time, event=event, X=X, S=S)
  dffull <- dffull[order(dffull[,"time"],decreasing=FALSE),]
  
  list_tbs <- list()
  unique_strat <- sort(unique(S))
  for(k in 1:length(unique_strat)){
    
    df <- dffull %>% filter(S==unique_strat[k])
    get_group <- df %>% group_by(across(starts_with("X"))) %>% summarise(count=n())
    cov_group <- get_group[,-ncol(get_group)]
    
    rtime_table1 <- evtime_table1 <- NULL
    intervals <- change.points[[k]]
    for(j in 2:length(intervals)){
      
      # table for risk times
      int <- intervals[j] - intervals[j-1]
      df$r <- ifelse(df$time > intervals[j], int, df$time-intervals[j-1])
      df$r[df$r < 0] <- 0
      tab <- df %>% group_by(across(starts_with("X"))) %>% summarise(risk_time=sum(r))
      rtime_table1 <- cbind(rtime_table1, tab$risk_time)
      
      #df$ev <- ifelse(df$time>=intervals[j-1] & df$time<intervals[j] & df$event==1, 1, 0)
      df$ev <- ifelse(df$time>intervals[j-1] & df$time<=intervals[j] & df$event==1, 1, 0)
      tab <- df %>% group_by(across(starts_with("X"))) %>% summarise(events_sum=sum(ev))
      evtime_table1 <- cbind(evtime_table1, tab$events_sum)
    }
    colnames(rtime_table1) <- 1:(ncol(rtime_table1))
    colnames(evtime_table1) <- 1:(ncol(evtime_table1))
    
    
    
    # table for events
    
    rtime_table <- cbind(cov_group, as.data.frame(rtime_table1))
    rtime_long <- gather(rtime_table, interval, rtime, (P+1):(P+ncol(rtime_table1)))
    evtime_table <- cbind(cov_group, as.data.frame(evtime_table1))
    evtime_long <- gather(evtime_table, interval, evtime, (P+1):(P+ncol(evtime_table1)))
    rtime_long$evt <- evtime_long$evtime
    table <- rtime_long %>% filter(rtime !=0)
    table$interval <- as.integer(table$interval)
    
    list_tbs[[k]] <- as.matrix(table)
    
  }
  
  return(list_tbs)
}

# function for creating current and historical tables
PCH <- function(time, event, X, S, historical, n.intervals, change.points, dCurrent){
  
  # create a list of intervals for the strata
  #change.points <- create_intervals(time, event, S, historical, n.intervals)
  
  # create tables for current data 
  if(dCurrent==TRUE){
    curr_tables <- suppressMessages(create_tables_new(time, event, X, S, change.points))
  }
  
  
  # create tables for historical data
  hist_tables <- list()
  if(length(historical)!=0){
    for(i in 1:length(historical)){
      data <- historical[[i]]
      tabs <- suppressMessages(create_tables_new(data[["time"]],data[["event"]],data[["X"]],S=data[["S"]],change.points))
      hist_tables <- c(hist_tables, list(tabs))
    }
  }
  
  if(dCurrent==FALSE){ 
    result <- list(hist_tables=hist_tables)
  }else if(length(historical)==0){ 
    result <- list(curr_tables=curr_tables) 
  }else{
    result <- list(curr_tables=curr_tables, hist_tables=hist_tables)
  }
  
  return(result)
  
}
  
  







#' Power/type I error calculation for the proportional hazards model with piecewise constant hazard and fixed a0
#'
#' @description Power/type I error calculation using power priors for the proportional hazards model with piecewise constant hazard and fixed \eqn{a_0}
#'
#' @param historical List of historical dataset(s). East historical dataset is stored in a list which contains four \emph{named} elements: \code{time}, \code{event}, \code{X} and \code{S}.
#' \itemize{
#' \item \code{time} is a vector of follow up times.
#' \item \code{event} is a vector of status indicators. Normally 0=alive and 1=dead.
#' \item \code{X} is a matrix of covariates. The first column must be the treatment indicator. 
#' \item \code{S} is a vector of integers, where each integer represents the stratum that the subject belongs to. For example, if there are three strata, S can take values 1, 2 or 3. 
#' }
#' @param a0 Vector containing numbers between 0 and 1 indicating the discounting parameter value for each historical dataset. The length of the vector should be equal to the length of \code{historical}.
#' @param n.subjects Number of subjects enrolled. 
#' @param n.events Number of events at which the trial will stop. 
#' @param n.intervals Vector of integers, indicating the number of intervals for the baseline hazards for each stratum. The length of the vector should be equal to the total number of strata. 
#' @param change.points List of vectors. Each vector in the list contains the change points for the baseline hazards for each stratum. The length of the vector should be equal to the total number of strata. 
#' By default, we assign the change points so that the same number of events are observed in all the intervals in the historical data.  
#' @param shared.blh Logical value indicating whether baseline hazard parameters are shared between the current and historical data. If TRUE, baseline hazard parameters are shared. 
#' @param samp.prior.beta Matrix of possible values of \eqn{\beta} to sample (with replacement) from. Each row is a possible \eqn{\beta} vector (a realization from the sampling prior for \eqn{\beta}).
#' @param samp.prior.lambda List of matrices, where each matrix represents the sampling prior for the baseline hazards for each stratum. The number of columns of each matrix should be equal to the number of intervals for that stratum.
#' @param x.samples (Only applies when there is no historical dataset) matrix of possible values of covariates from which covariate vectors are sampled with replacement.   
#' @param s.samples (Only applies when there is no historical dataset) vector of possible values of the stratum index from which the stratum indices are sampled with replacement.   
#' @param dist.enroll Distribution for enrollment times. The choices are "Uniform" or "Exponential". 
#' @param param.enroll Parameter for the distribution of enrollment times. If \code{dist.enroll} is "Uniform", the enrollment times follow Unif(0, \code{param.enroll}). If \code{dist.enroll} is "Exponential", 
#' the enrollment times follow Exponential(rate=\code{param.enroll}).
#' @param rand.prob Randomization probability for the treated group. The default value is 0.5.  
#' @param prob.drop Probability of subjects dropping out of the study (non-administrative censoring). The default value is zero.  
#' @param param.drop Parameter for dropout time simulations. The dropout times follow Unif(0, \code{param.drop}). The default value is zero. 
#' @param dist.csr Distribution for (administrative) censorship times. The choices are "Uniform", "Constant" and "Exponential". The default choice is "Constant". 
#' @param param.csr Parameter for the (administrative) censorship times. If \code{dist.csr} is "Uniform", the censorship times follow Unif(0, \code{param.csr}). 
#' If \code{dist.csr} is "Constant", the censorship times of all subjects are equal to \code{param.csr}.
#' If \code{dist.csr} is "Exponential", the censorship times follow Exponential(rate=\code{param.csr}).
#' The default value is 10^4.  
#' @param min.follow.up Minimum amount of time for which subjects are followed up. The default value is zero. 
#' @param max.follow.up Maximum amount of time for which subjects are followed up. The default value is 10^4.  
#' @param prior.beta Prior used for \eqn{\beta}. The choices are "Uniform" and "Normal". If \code{prior.beta} is "Uniform", the uniform improper prior is used.
#' If \code{prior.beta} is "Normal", independent normal priors are used for each element of \eqn{\beta}. The default choice is "Normal". 
#' @param prior.beta.mean (Only applies if \code{prior.beta} is "Normal") vector of means of the normal prior on \eqn{\beta}. The default value is zero for all the elements of \eqn{\beta}. 
#' @param prior.beta.sd (Only applies if \code{prior.beta} is "Normal") vector of standard deviations of the normal prior on \eqn{\beta}. The default value is 10^3 for all the elements of \eqn{\beta}. 
#' @param prior.lambda Prior used for \eqn{\lambda}. The choices are "Gamma", "Log-normal" and "Improper". The default choice is "Gamma". 
#' 
#' If \code{prior.lambda} is "Gamma", then the prior on the first element of \eqn{\lambda} is 
#' 
#' Gamma(shape=\code{prior.lambda.hp1[1]}, rate=\code{prior.lambda.hp2[1]}).
#' 
#' If \code{prior.lambda} is "Log-normal", then the prior on the first element of \eqn{\lambda} is Log-normal(mean=\code{prior.lambda.hp1[1]}, sd=\code{prior.lambda.hp2[1]}).
#' 
#' If \code{prior.lambda} is "Improper", then the prior on each element of \eqn{\lambda} is the improper prior \eqn{\lambda^{-1}}. 
#' @param prior.lambda.hp1 (Only applies if \code{prior.lambda} is "Gamma" or "Log-normal") Vector of first hyperparameters of the prior on \eqn{\lambda}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda}.
#' @param prior.lambda.hp2 (Only applies if \code{prior.lambda} is "Gamma" or "Log-normal") Vector of second hyperparameters of the prior on \eqn{\lambda}. 
#' The length of the vector should be equal to the dimension of \eqn{\lambda}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda}.
#' @param prior.lambda0.hp1 (Only applies if \code{shared.blh} is FALSE and if \code{prior.lambda} is "Gamma" or "Log-normal") Vector of first hyperparameters of the prior on \eqn{\lambda_0}.
#' We assume the same distribution choice for the prior for \eqn{\lambda_0} and \eqn{\lambda}.
#' The length of the vector should be equal to the dimension of \eqn{\lambda_0}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda_0}.
#' @param prior.lambda0.hp2 (Only applies if \code{shared.blh} is FALSE and if \code{prior.lambda} is "Gamma" or "Log-normal") Vector of second hyperparameters of the prior on \eqn{\lambda_0}.  
#' We assume the same distribution choice for the prior for \eqn{\lambda_0} and \eqn{\lambda}.
#' The length of the vector should be equal to the dimension of \eqn{\lambda_0}, i.e., the total number of intervals for all strata. The default value is 10^(-5) for all the elements of \eqn{\lambda_0}. 
#' @param lower.limits Vector of lower limits for parameters (\eqn{\beta}, \eqn{\lambda}, and \eqn{\lambda_0}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is -100 for \eqn{\beta} and 0 for \eqn{\lambda} (may not be appropriate for all situations). 
#' @param upper.limits Vector of upper limits for parameters (\eqn{\beta}, \eqn{\lambda}, and \eqn{\lambda_0}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is 100 for all parameters (may not be appropriate for all situations).
#' @param slice.widths Vector of initial slice widths for parameters (\eqn{\beta}, \eqn{\lambda}, and \eqn{\lambda_0}, in this order) to be used by the slice sampler. The length of the vector should be equal to the total number of parameters. The default is 0.1 for all parameters (may not be appropriate for all situations). 
#' @param nMC Number of iterations (excluding burn-in samples) for the slice sampler. The default is 10,000.
#' @param nBI Number of burn-in samples for the slice sampler. The default is 250.
#' @param delta Prespecified constant that defines the boundary of the null hypothesis. The default is zero.
#' @param nullspace.ineq Character string specifying the inequality of the null hypothesis. The options are ">" and "<". If ">" is specified, the null hypothesis is \eqn{H_0}: \eqn{\beta_1} \eqn{\ge} \eqn{\delta}. If "<" is specified, the null hypothesis is \eqn{H_0}: \eqn{\beta_1} \eqn{\le} \eqn{\delta}. The default choice is ">".
#' @param gamma Posterior probability threshold for rejecting the null. The null hypothesis is rejected if posterior probability is greater \code{gamma}. The default is 0.95.
#' @param N Number of simulated datasets to generate. The default is 10,000.
#'
#' @details The proportional hazards model with piecewise constant hazard is implemented. 
#' We assume \eqn{\beta} is the regression coefficients. We assume the first column of the covariate matrix is the treatment indicator,
#' and the corresponding parameter is \eqn{\beta_1}. The baseline hazards of the current data are denoted by \eqn{\lambda}. 
#' If the baseline hazards are not shared between the historical and current data, 
#' the baseline hazards for the historical data are denoted by \eqn{\lambda_0}.
#' 
#' To perform sample size determination, we test the hypotheses
#' \deqn{H_0: \beta_1 \ge \delta} and \deqn{H_1: \beta_1 < \delta.}
#' 
#' If historical datasets are provided, the algorithm samples with replacement from the historical covariates to construct the simulated datasets.
#' Otherwise, the algorithm samples with replacement from \code{x.samples}. One of the arguments \code{historical} and \code{x.samples} must be provided.
#'
#' The sampling prior for the treatment parameter can be generated from a normal distribution (see examples). 
#' For example, suppose one wants to compute the power for the hypotheses \eqn{H_0: \beta_1 \ge 0} and \eqn{H_1: \beta_1 < 0.} 
#' To approximate the sampling prior for \eqn{\beta_1}, one can simply sample from a normal distribution with negative mean, 
#' so that the mass of the prior falls in the alternative space. Conversely, to compute the type I error rate, one can 
#' sample from a normal distribution with positive mean, so that the mass of the prior falls in the null space. 
#' 
#' Posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{\lambda}. The default upper limits 
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' If a sampling prior with support in the null space is used, the value returned is a Bayesian type I error rate.
#' If a sampling prior with support in the alternative space is used, the value returned is a Bayesian power.
#'
#' 
#' @return Power or type I error is returned, depending on the sampling prior used.
#' The posterior probabilities of the alternative hypothesis are returned.  
#' The average posterior mean of \eqn{\beta} is also returned. 
#' @references Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{phm.fixed.a0}}
#' @examples
#' 
#' # Simulate two historical datasets
#' set.seed(1)
#' n <- 100
#' P <- 4
#' time1 <- round(rexp(n, rate=0.5),1)
#' event1 <- rep(1,n)
#' X1 <- matrix(rbinom(n*P,p=0.5,size=1), ncol=P)
#' S1 <- c(rep(1,n/2),rep(2,n/2))
#' time2 <- round(rexp(n, rate=0.7),1)
#' event2 <- rep(1,n)
#' X2 <- matrix(rbinom(n*P,p=0.5,size=1), ncol=P)
#' S2 <- c(rep(1,n/2),rep(2,n/2))
#' historical <- list(list(time=time1, event=event1, X=X1, S=S1),
#'                    list(time=time2, event=event2, X=X2, S=S2))
#' 
#' a0 <- c(0.3, 0.6)
#' n.subjects <- 100
#' n.events <- 30
#' 
#' # We choose three intervals for the first stratum and two intervals for the second stratum
#' n.intervals <- c(3,2) 
#' change.points <- list(c(1,2), 2)
#' 
#' # Generate sampling priors
#' 
#' # The null hypothesis here is H0: beta_1 >= 0. To calculate power,
#' # we can provide samples of beta_1 such that the mass of beta_1 < 0.
#' # To calculate type I error, we can provide samples of beta_1 such that
#' # the mass of beta_1 >= 0.
#' samp.prior.beta1 <- rnorm(100, mean=-1, sd=1)
#' # Here, mass is put on the alternative region, so power is calculated.
#' samp.prior.beta <- cbind(samp.prior.beta1, matrix(rnorm(100*(P-1)), 100, P-1))
#' 
#' # Point mass sampling priors are used for lambda
#' lambda_strat1 <- matrix(c(0.5, 0.5, 0.5), nrow=1)
#' lambda_strat2 <- matrix(c(0.7, 0.7), nrow=1)
#' samp.prior.lambda <- list(lambda_strat1, lambda_strat2)
#' 
#' 
#' nMC <- 100 # nMC should be larger in practice
#' nBI <- 50
#' N <- 5 # N should be larger in practice
#' 
#' 
#' result <- power.phm.fixed.a0(historical=historical, a0=a0, n.subjects=n.subjects, 
#'                              n.events=n.events, n.intervals=n.intervals,  
#'                              change.points=change.points, shared.blh=FALSE,
#'                              samp.prior.beta=samp.prior.beta, 
#'                              samp.prior.lambda=samp.prior.lambda,
#'                              dist.enroll="Uniform", param.enroll=0.5,
#'                              nMC=nMC, nBI=nBI, delta=0, nullspace.ineq=">", N=N)
#' result$`power/type I error`
#' 
#'
#' @export
power.phm.fixed.a0 <- function(historical, a0, n.subjects, n.events, 
                               n.intervals, change.points=NULL, shared.blh,
                               samp.prior.beta, samp.prior.lambda, # list of matrices
                               x.samples=matrix(), s.samples=NULL, # these two are matrices
                               dist.enroll, param.enroll,
                               rand.prob=0.5, prob.drop=0, param.drop=0, 
                               dist.csr="Constant", param.csr=10000, min.follow.up=0, max.follow.up=10000,
                               prior.beta="Normal", prior.beta.mean=rep(0,50), prior.beta.sd=rep(1000,50),
                               prior.lambda="Gamma", prior.lambda.hp1=rep(10^(-5),50), prior.lambda.hp2=rep(10^(-5),50), 
                               prior.lambda0.hp1=rep(10^(-5),50), prior.lambda0.hp2=rep(10^(-5),50),
                               lower.limits=NULL, upper.limits=rep(100, 50), slice.widths=rep(0.1, 50),
                               nMC=10000, nBI=250, delta=0, nullspace.ineq=">", gamma=0.95, N=10000){
  
  # add zero and infinity to change.points
  change.points.new <- list()
  if(is.null(change.points)){
    change.points.new <- create_intervals_design(historical,n.intervals)
  }else{
    for(i in 1:length(change.points)){
      l <- change.points[[i]]
      l1 <- unique(c(0, l, Inf))
      change.points.new[[i]] <- l1
    }
  }


  # build matrix of covariates to sample from
  P <- ncol(historical[[1]]$X)
  if(length(historical)==0){
    x_sim <- x.samples; # dimension P-1
  }else{
    x_sim <- matrix(NA, nrow=0, ncol=P-1)
    for(k in 1:length(historical)){
      dat <- historical[[k]]
      x_h <- as.matrix(dat$X[,-1])
      x_sim <- rbind(x_sim, x_h)
    }
  }
  # build matrix of strata to sample from
  if(length(historical)==0){
    s_sim <- s.samples; # dimension 1
  }else{
    s_sim <- NULL
    for(k in 1:length(historical)){
      dat <- historical[[k]]
      s_h <- dat$S
      s_sim <- c(s_sim, s_h)
    }
  }
  
  # repeat N times 
  beta_mean <- matrix(NA, nrow=N, ncol=P)
  #lambda_mean1 <- NULL
  #lambda_mean2 <- NULL
  posterior_probs <- NULL
  for(l in 1:N){
    # print(l)
    # sample stratum variable
    ind <- sample(1:length(s_sim), n.subjects, replace = TRUE)
    s <- s_sim[ind]
    
    # sample beta
    ind <-  sample(1:nrow(samp.prior.beta), 1, replace = TRUE)
    beta <- samp.prior.beta[ind,]
    
    # sample lambda
    list_lambda <- list()
    for(j in 1:length(unique(s))){
      lams <- samp.prior.lambda[[j]]
      ind <-  sample(1:nrow(lams), 1, replace = TRUE)
      list_lambda <- c(list_lambda, list(lams[ind,]))
    }
    
    # sample enrollment times
    if(dist.enroll=="Uniform"){
      r <- runif(n.subjects,0,param.enroll)
    }else{
      r <- rexp(n.subjects,rate=param.enroll)
    }
    
    # sample covariates and treatment variable
    ind <-  sample(1:nrow(x_sim), n.subjects, replace = TRUE)
    x <- x_sim[ind,]
    z <- rbinom(n.subjects,size=1,rand.prob)
    x <- cbind(trt=z, x)
    
    
    phi <- exp(x%*%beta)
    
    # simulate event time ti
    eventtimes <- NULL

    for(i in 1:n.subjects){

      strat <- s[i]
      lambda_s <- list_lambda[[strat]]
      st <- change.points.new[[strat]]
      st_1 <- st[-length(st)]

      
      k <- 1
      
      while(TRUE){

        theta_k <- phi[i]*lambda_s[k]
        
        t_tilde <- rexp(1,rate=theta_k) + st_1[k]
        if(k > length(st_1) | t_tilde <= st[k+1]){
          break
        }
        k <- k + 1
      }
      eventtimes <- c(eventtimes, t_tilde)
    }
    
    # sample (administrative) censorship times
    if(dist.csr=="Uniform"){
      ctime <- runif(n.subjects,0,param.csr)
    }else if(dist.csr=="Constant"){
      ctime <- rep(param.csr, n.subjects)
    }else{
      ctime <- rexp(n.subjects,rate=param.csr)
    }
    y <- ifelse(eventtimes <= ctime, eventtimes, ctime)
    nu <- ifelse(eventtimes <= ctime, 1, 0)
    
    # simulate dropouts
    # simulate a dropout time; if dropout time < y, the person is a dropout; 
    # otherwise the person has an event or is administratively censored
    num_drops <- round(n.subjects * prob.drop, 0)
    
    if(num_drops > 0){
      drop_ind <- sample(1:n.subjects, size=num_drops, replace=FALSE)
      droptime <- runif(num_drops, 0, param.drop)
      obstime <- ifelse(y[drop_ind] <= droptime, y[drop_ind], droptime)
      bool_drop <- ifelse(y[drop_ind] <= droptime, 0, 1)
      nu[drop_ind][bool_drop==1] <- 0
      y[drop_ind] <- obstime
    }
    
    e <- r+y
    complete_data <- data.frame(X=x, S=s, enrtimes=r, eventtimes=eventtimes, ctime=ctime, y=y, nu=nu, t_elps=e)

    # create original data
    df_events <- complete_data[complete_data$nu==1,]
    df_events1 <- df_events[order(df_events$t_elps),]
    if(nrow(df_events1) < n.events){
      stoptime <- max.follow.up
    }else{
      stoptime <- df_events1[n.events, "t_elps"]
      if(stoptime < min.follow.up){
        stoptime <- min.follow.up
      }
    }

    finaldf <- complete_data[complete_data$enrtimes < stoptime,]
    finaldf$new_y <- ifelse(finaldf$t_elps > stoptime, stoptime-finaldf$enrtimes, finaldf$y)
    finaldf$new_nu <- ifelse(finaldf$t_elps > stoptime, 0, finaldf$nu)
    
    # create tables
    tables <- PCH(time=finaldf$new_y, event=finaldf$new_nu, X=finaldf[,1:ncol(x)], S=finaldf$S, 
                  historical=historical, n.intervals=n.intervals, change.points=change.points.new, dCurrent=TRUE)
    t1 <- tables[["curr_tables"]]
    t2 <- tables[["hist_tables"]]
    #print(t1)
    #print(t2)
    if(is.null(lower.limits)){
      lower.limits = c(rep(-100,P), rep(0,50))
    }
    samples <- phm_fixed_a0(t1, t2, a0,
                            n.intervals, shared.blh, P, prior.beta, prior.beta.mean, prior.beta.sd,
                            prior.lambda, prior.lambda.hp1, prior.lambda.hp2, 
                            prior.lambda0.hp1, prior.lambda0.hp2,
                            lower.limits, upper.limits, slice.widths, 
                            nMC, nBI, dCurrent=TRUE)
    
    # compute probability of success
    beta1 <- samples$beta_samples[,1]
    
    if(nullspace.ineq == ">"){
      pr <- sum(beta1 < delta) / length(beta1)
    }else{
      pr <- sum(beta1 > delta) / length(beta1)
    }
    
    posterior_probs <- c(posterior_probs, pr)
    
    beta_mean[l,] <- colMeans(samples$beta_samples)
    
    
    #lambda_mean1 <- c(lambda_mean1, colMeans(samples$lambda_samples[[1]]))
    #lambda_mean2 <- c(lambda_mean2, colMeans(samples$lambda_samples[[2]]))
    
    #finaldf$dataset_id <- l
    #simdf <- rbind(simdf, finaldf)
    
  }
  
  power <- mean(posterior_probs >= gamma)
  return(list(#"simulated dataset"=simdf, 
              "power/type I error"=power, 
              "posterior probabilities"=posterior_probs,
              "average posterior mean of beta"=beta_mean))
              #"lambda_mean1"=lambda_mean1,"lambda_mean2"=lambda_mean2))
}





#' Model fitting for the proportional hazards model with piecewise constant hazard and fixed a0
#'
#' @description Model fitting using power priors for the proportional hazards model with piecewise constant hazard and fixed a0
#' @param time Vector of follow up times.
#' @param event Vector of status indicators. Normally 0=alive and 1=dead.
#' @param X Matrix of covariates. The first column must be the treatment indicator. 
#' @param S Vector of integers, where each integer represents the stratum that the subject belongs to. For example, if there are three strata, S can take values 1, 2 or 3. 
#' @param change.points List of vectors. Each vector in the list contains the change points for the baseline hazards for each stratum. The length of the vector should be equal to the total number of strata. 
#' By default, we assign the change points so that the same number of events are observed in all the intervals in the pooled current and historical data.  
#' @inheritParams power.phm.fixed.a0
#' 
#' @details The proportional hazards model with piecewise constant hazard is implemented. 
#' We assume \eqn{\beta} is the regression coefficients. We assume the first column of the covariate matrix is the treatment indicator,
#' and the corresponding parameter is \eqn{\beta_1}. The baseline hazards of the current data are denoted by \eqn{\lambda}. 
#' If the baseline hazards are not shared between the historical and current data, 
#' the baseline hazards for the historical data are denoted by \eqn{\lambda_0}.
#' 
#' Posterior samples are obtained through slice sampling. 
#' The default lower limits are -100 for \eqn{\beta} and 0 for \eqn{\lambda}. The default upper limits 
#' for the parameters are 100. The default slice widths for the parameters are 0.1.
#' The defaults may not be appropriate for all situations, and the user can specify the appropriate limits
#' and slice width for each parameter.
#'
#' 
#' @return Posterior samples of \eqn{\beta}, \eqn{\lambda} and \eqn{\lambda_0} (if baseline hazards are not shared between the current and historical data) are returned.
#' @references Neal, Radford M. Slice sampling. Ann. Statist. 31 (2003), no. 3, 705--767.
#' @seealso \code{\link{power.phm.fixed.a0}}
#' @examples 
#' 
#' set.seed(1)
#' # Simulate current data
#' n <- 50
#' P <- 4
#' time <- round(rexp(n, rate=0.5),1)
#' event <- rep(1,n)
#' X <- matrix(rbinom(n*P,p=0.5,size=1), ncol=P)
#' S <- c(rep(1,n/2),rep(2,n/2))
#' 
#' # Simulate two historical datasets
#' n <- 100
#' time1 <- round(rexp(n, rate=0.5),1)
#' event1 <- rep(1,n)
#' X1 <- matrix(rbinom(n*P,p=0.5,size=1), ncol=P)
#' S1 <- c(rep(1,n/2),rep(2,n/2))
#' time2 <- round(rexp(n, rate=0.7),1)
#' event2 <- rep(1,n)
#' X2 <- matrix(rbinom(n*P,p=0.5,size=1), ncol=P)
#' S2 <- c(rep(1,n/2),rep(2,n/2))
#' historical <- list(list(time=time1, event=event1, X=X1, S=S1),
#'                    list(time=time2, event=event2, X=X2, S=S2))
#' a0 <- c(0.3, 0.6)
#' 
#' 
#' # We choose three intervals for the first stratum and two intervals for the second stratum
#' n.intervals <- c(3,2) 
#' change.points <- list(c(1,2), 2)
#' 
#' 
#' nMC <- 1000 # nMC should be larger in practice
#' nBI <- 50
#' 
#' result <- phm.fixed.a0(time=time, event=event, X=X, S=S,
#'                        historical=historical, a0=a0, n.intervals=n.intervals, 
#'                        change.points=change.points, shared.blh=FALSE, 
#'                        nMC=nMC, nBI=nBI)
#' 
#' colMeans(result$beta_samples) # posterior mean of beta
#' colMeans(result$lambda_samples[[1]]) # posterior mean of baseline hazards for stratum 1
#' colMeans(result$lambda_samples[[2]]) # posterior mean of baseline hazards for stratum 2
#'
#' @export
phm.fixed.a0 <- function(time, event, X, S, historical, a0, n.intervals, change.points=NULL, shared.blh,
                        prior.beta="Normal", prior.beta.mean=rep(0,50), prior.beta.sd=rep(1000,50),
                        prior.lambda="Gamma", prior.lambda.hp1=rep(10^(-5),50), prior.lambda.hp2=rep(10^(-5),50), 
                        prior.lambda0.hp1=rep(10^(-5),50), prior.lambda0.hp2=rep(10^(-5),50),
                        lower.limits=NULL, upper.limits=rep(100, 50), slice.widths=rep(0.1, 50),
                        nMC=10000, nBI=250){
  
  # add zero and infinity to change.points
  change.points.new <- list()
  if(is.null(change.points)){
    change.points.new <- create_intervals(time, event, S, historical, n.intervals)
  }else{
    for(i in 1:length(change.points)){
      l <- change.points[[i]]
      l1 <- unique(c(0, l, Inf))
      change.points.new[[i]] <- l1
    }
  }
  # create tables
  tables <- PCH(time=time, event=event, X=X, S=S, 
                historical=historical, n.intervals=n.intervals, change.points=change.points.new, dCurrent=TRUE)
  t1 <- tables[["curr_tables"]]
  t2 <- tables[["hist_tables"]]

  if(is.null(lower.limits)){
    lower.limits = c(rep(-100,P), rep(0,50))
  }
  
  P <- ncol(X)
  samples <- phm_fixed_a0(t1, t2, a0,
                          n.intervals, shared.blh, P, prior.beta, prior.beta.mean, prior.beta.sd,
                          prior.lambda, prior.lambda.hp1, prior.lambda.hp2, 
                          prior.lambda0.hp1, prior.lambda0.hp2,
                          lower.limits, upper.limits, slice.widths, 
                          nMC, nBI, dCurrent=TRUE)
  return(samples)
}
    
