# FILE WITH ALL THE NEEDED PROGRAMS FOR ESTIMATION OF R0 #
#  ALSO CONTAINS FUNCTIONS FOR SIMULATING DATA           #
##########################################################
# CONTENTS #
# 1.  getR0SI.WP(N,distn="gamma",k=6)
# 2.  getR0.WP(N,pars,distn="gamma",k=20)
# 3.  getRt.WP(Ns,ps)
# 4.  getRt.WT(Ns,ps)
# 5.  likN.MN(params,N)
# 6.  likN.gam(params,N,k)
# 7.  likN.logit.SIknown(params,distn="gamma",dist.pars,N,k=6)
# 8.  likN.logit.SIunknown(params,N,k=6)
# 9.  getNsMN(R0,ps,num.samples=1000,epi.size,num.start,perc.off=0.02)
# 10.  getNsGam(R0,shape,rate,num.samples=1000,epi.size,num.start)
# 11. getNsMNL(R0,ps,num.samples=500,max.days,num.start)
# 12. getNsGamL(R0,shape,rate,num.samples=500,max.days,num.start)
# 13. getSummaryMN(pars): summary stats for MN Serial interval
# 14. covProbs(parsEsts,k,epi.size,num.start,numRep=1000,distn="MN",est.si=T): coverage probabilities for WP method estimates

##########################################
# 1. FUNCTION TO GET R0 WHEN SI IS UNKNOWN #
##########################################
getR0SI.WP <- function(N,distn="gamma",k=6){

	if(distn=="gamma"){
		tmp <- constrOptim(c(1.2,4,0.5),likN.gam,N=N,k=k,ui=diag(1,3),ci=rep(0,3),method="Nelder-Mead")
		params.tmp <- tmp$par
		
		if(tmp$convergence>0){return(c("ERROR: did not converge"))}
		if(tmp$convergence==0){return(params.tmp)}
		
	}		

	if(distn=="MN"){
		tmp <- constrOptim(rep(1.2,k),likN.MN,N=N,ui=diag(1,k),ci=rep(0,k),method="Nelder-Mead")
		params.tmp <- tmp$par
		ps  <- params.tmp/sum(params.tmp)
		R0 <- sum(params.tmp)
		if(tmp$convergence>0){return(c("ERROR: did not converge"))}
		if(tmp$convergence==0){return(c(R0,ps))}
	}

}

##########################################
# 2. FUNCTION TO GET R0 WHEN SI IS KNOWN #
##########################################
getR0.WP <- function(N,pars,distn="gamma",k=20){
# function to get R0 when the serial interval is known
    ps <- rep(0,k)

    # get the si, depending on how it is distributed
    if(distn=="gamma"){
        alpha <- pars[1]
        rate <- pars[2]
        for(i in 1:k){
            ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
        }
    }
    if(distn=="exp"){
        for(i in 1:k){
            ps[i] <- pexp(i,pars)-pexp(i-1,pars)
        }
    }
    if(distn=="MN"){
        ps <- pars
    }
    
    ps <- ps/sum(ps)
    
    N <- N[1:max(which(N!=0))] # cut off zero values at the end (since some values at the beginning might be zeroed out)
    
    N.long <- c(rep(0,(k-1)),N)

    l <- length(N.long)

    if(l > k){    
        sumK <- 0
        for(i in k:(l-1)){
            sumK <- ps%*%N.long[i:(i-k+1)] + sumK
        }
        R0 <- sum(N)/sumK
    } else R0 <- NA    

    return(R0)
}

####################################################
# PROGRAMS TO GET ESTIMATES OF R_T                 #
####################################################
# 3 #
getRt.WP <- function(N,dist.pars,k=6){

    if(distn=="gamma"){
        alpha <- dist.pars[1]
        rate <- dist.pars[2]
        for(i in 1:k){
            ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
        }
    
    ps <- ps/sum(ps)
    }
    if(distn=="MN"){
        ps <- dist.pars
    }

	ests <- constrOptim(rep(1.2,k),likN.logit.SIknown,distn="MN",dist.pars=ps,N=N,ui=diag(1,k),ci=rep(0,k),method="Nelder-Mead",k=k)	
	params <- ests$par

	times <- seq(1,length(N),1)

	getRt <- function(pars,times){return(pars[1]+pars[2]/(1+exp(pars[3]*(times-pars[4]))))}

	Rts <- getRt(params,times)

	return(Rts)
}

# 4. wallinga and teunis #
getRt.WT <- function(N,distn="gamma",dist.pars,k=6){
	# ps are the estimate of the serial interval
	ps <- rep(0,k)
    if(distn=="gamma"){
        alpha <- dist.pars[1]
        rate <- dist.pars[2]
        for(i in 1:k){
            ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
        }
    
    ps <- ps/sum(ps)
    }
    if(distn=="MN"){
        ps <- dist.pars
    }

	LenN <- length(N)

    	N <- c(rep(0,k),N,rep(0,k))
	Rs <- rep(0,LenN) 

	for(j in (k+1):(LenN+k)){
		if(N[j] >0){

		tmp <- 0
		for (jj in 1:k){
            	denom <- N[(j+jj-1):(j+jj-k)]%*%ps
            	if(denom != 0 & !is.na(denom)){tmp <- N[j+jj]*ps[jj]/denom + tmp}
			else {tmp <- N[j+jj]}
		}
    
		Rs[j-k] <- tmp
		}else Rs[j-k] <- NA
	}
	
	return(Rs)

}


########################
# LIKELIHOOD FUNCTIONS #
########################

# N is the vector of daily case counts
# k is the maximal length of the serial interval
# alpha and rate are the shape and rate of the gamma distrn
# in all the programs, thetas=R0*ps
# lambda is sum_j p_j N_{t-j}

############################################################
# LIKELIHOODS FOR EXPONENTIAL GROWTH PHASE OF THE EPIDEMIC #
############################################################
# 5. Likelihood for multinomial SI
#####################
likN.MN <- function(params,N){
#params is of the form (thetas), where theta_i=R0 p_i

    k <- length(params)
    thetas <- params   

    N.long <- c(rep(0,(k-1)),N)
    l <- length(N.long)
        
    lik <- 0
    
    for(i in k:(l-1)){
        lambda <- thetas%*%N.long[i:(i-k+1)] 
        lik <- -lambda + N.long[i+1]*log(lambda) + lik
    }
    
    return(-lik)
}

###################
# 6. Likelihood for gamma SI
#####################
likN.gam <- function(params,N,k){
#params is of the form (R0,shape,rate)

    R0 <- params[1]
    alpha <- params[2]
    rate <- params[3]

    ps <- rep(0,k)
    
    for(i in 1:k){
        ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
    }
    ps <- ps/sum(ps)
    
    thetas <- R0*ps

    N.long <- c(rep(0,(k-1)),N)
    l <- length(N.long)
        
    lik <- 0
    
    for(i in k:(l-1)){
        lambda <- thetas%*%N.long[i:(i-k+1)] 
        lik <- -lambda + N.long[i+1]*log(lambda) + lik
    }
    
    return(-lik)
}

#######################################################
# LIKELIHOODS FOR THE ENTIRE EPIDEMIC,                #
# ASSUMING THAT R_t FOLLOWS A FOUR PARAMETER LOGISTIC #
#######################################################
## Likelihoods to estimate R_t throughout the crouse of an epidemic.
## We assume that R_t follows a four parameter logistic.
## Serial interval is a multinomial distribution.

# 7
likN.logit.SIknown <- function(params,distn="gamma",dist.pars,N,k=6){
# likelihood for the model that uses a four parameter logistic to fit the change in R through time
#  params is of the form c(a,b,c,d)
#  dist.pars=(alpha,rate) if distn="gamma" or (ps) is distn="MN"
#  We assume that the serial is known in this case

    a <- params[1]
    b <- params[2]
    c <- params[3]
    d <- params[4]
    
    if(distn=="gamma"){
        alpha <- dist.pars[1]
        rate <- dist.pars[2]
        for(i in 1:k){
            ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
        }
    
    ps <- ps/sum(ps)
    }
    if(distn=="MN"){
        ps <- dist.pars
    }
    
  
    N.long <- c(rep(0,(k-1)),N)
    l <- length(N.long)
        
    lik <- 0

    for(i in k:(l-1)){    
        pN <- ps%*%N.long[i:(i-k+1)] # ps*previous Ns
        R.t <- a+b/(1+exp(c*(i-k+1-d))) # value of R_t for i, we assume that i=k corresponds to t=0

        lambda <- R.t*pN
            
        lik <- -lambda + N.long[i+1]*log(lambda) + lik #the value of the likelihood
    }
    
    return(-lik)
}

# 8 #
likN.logit.SIunknown <- function(params,N,k=6){
# likelihood for the model that uses a four parameter logistic to fit the change in R through time
#  params is of the form c(a,b,c,d,alpha,beta)
#  We assume that the serial is unknown in this case
#  the serial interval must be gamma distributed (I can't figure out how to estimate it with a MN)

    a <- params[1]
    b <- params[2]
    c <- params[3]
    d <- params[4]
    alpha <- params[5]
    rate <- params[6]

    for(i in 1:k){
        ps[i] <- pgamma(i,shape=alpha,rate=rate)-pgamma(i-1,shape=alpha,rate=rate)
    } 
    ps <- ps/sum(ps)
   
    N.long <- c(rep(0,(k-1)),N)
    l <- length(N.long)
        
    lik <- 0

    for(i in k:(l-1)){    
        pN <- ps%*%N.long[i:(i-k+1)] # ps*previous Ns
        R.t <- a+b/(1+exp(c*(i-k+1-d))) # value of R_t for i, we assume that i=k corresponds to t=0

        lambda <- R.t*pN
            
        lik <- -lambda + N.long[i+1]*log(lambda) + lik #the value of the likelihood
    }
    
    return(-lik)
}

######################################
# PROGRAMS TO SIMULATE EPIDEMIC DATA #
######################################
# SIMULATE EPIS OF A GIVEN SIZE. LENGTH OF EPI IS NOT CONSTRAINED

# 9. MULTINOMIAL DISTRIBUTION
getNsMN <- function(R0,ps,num.samples=1000,epi.size,num.start,perc.off=0.02){
    # function to generate data for given parameter values with Multinomial serial interval. 
    # final epi size must be within perc.off of actual value.
    # R0, ps are the parameters from which we simulate.
    
set.seed(21)
library(Hmisc)
#library(Hmisc,lib.loc='~/Rlibs')

M <- epi.size # final size of the epidemic (number of cases)
bdd <- perc.off*M # amount the simulation can be off.

ps.mat <- t(as.matrix(ps))

# space allocation #
m <- 1
ind <- 1

while(m <= num.samples){
    day.counts <- new.cts.tab <- c(num.start)
    new.cases <- rep(1,num.start)
    new.days <- c(1)
    min.ind <- min(c(1:length(new.cts.tab))[new.cts.tab != 0]) 
    
    while(sum(day.counts[1:min.ind]) < (M+bdd)){
        num.new <- rpois(length(new.cts.tab),new.cts.tab*R0) # vector of the number of new cases generated by each day of cases (X_i dot)
    
        n <- sum(num.new) #total number of new cases on this iteration
        if(n == 0) break

        new.cases.start <- rep(new.days,num.new[new.days])
        new.cases <- new.cases.start + rMultinom(ps.mat,n) #cases that are new for this iteration
        new.cases <- sort(new.cases)[1:min(M,length(new.cases))] # get rid of unneeded cases (might be overly conservative)

        new.max <- max(new.cases,length(day.counts)) # maximal day in the data
        day.counts <- c(day.counts,rep(0,new.max-length(day.counts))) # 
        new.counts <- as.vector(table(new.cases))
        new.cts.tab <- rep(0,new.max)
        new.days <- unique(new.cases)
        new.cts.tab[sort(new.days)] <- new.counts
        
        day.counts <- day.counts+new.cts.tab  #allows to keep going if not close to the end

        min.ind <- min(which(new.cts.tab != 0))  
    }

    ind.max.dc <- max(which(cumsum(day.counts)<=(M+bdd)))
    sum.dc <- sum(day.counts[1:ind.max.dc])
    if((sum.dc <= (M+bdd)) & (sum.dc >= (M-bdd)))    day.counts <- day.counts[1:ind.max.dc]

    if((sum.dc <= (M+bdd)) & (sum.dc >= (M-bdd))){
        l.dc <- length(day.counts)
    
        if(m==1){
            Ns <- matrix(0,nr=1,nc=l.dc)
            Ns[m,] <- day.counts
        }
        if(m > 1){
            l.Ns <- length(Ns[m-1,])
            if(l.Ns > l.dc)  Ns <- rbind(Ns,c(day.counts,rep(0,l.Ns-l.dc)))
            if(l.Ns < l.dc)  Ns <- rbind(cbind(Ns,matrix(0,nc=(l.dc-l.Ns),nr=(m-1))),day.counts)
            if(l.Ns == l.dc) Ns <- rbind(Ns,day.counts)
        }
        
        m <- m+1
    }
    ind <- ind+1
 #   print(cbind(ind,m))
}
row.names(Ns) <- NULL

return(Ns)

}


# 10. GAMMA DISTRIBUTION
getNsGam <- function(R0,shape,rate,num.samples=1000,epi.size,num.start,perc.off=0.02){
    # function to generate data for given parameter values
    
set.seed(21)

M <- epi.size
bdd <- perc.off*M # amount the simulation can be off.

# space allocation #
m <- 1
ind <- 1

while(m <= num.samples){
    day.counts <- new.cts.tab <- c(num.start)
    new.cases <- rep(1,num.start)
    new.days <- c(1)
    min.ind <- min(c(1:length(new.cts.tab))[new.cts.tab != 0]) 
    
    while(sum(day.counts[1:min.ind]) < (M+bdd)){
        num.new <- rpois(length(new.cts.tab),new.cts.tab*R0) # vector of the number of new cases generated by each day of cases (X_i dot)
    
        n <- sum(num.new) #total number of new cases on this iteration
        if(n == 0) break

        new.cases.start <- rep(new.days,num.new[new.days])
        new.cases <- new.cases.start + ceiling(rgamma(n,shape=shape,rate=rate)) #cases that are new for this iteration
        new.cases <- sort(new.cases)[1:min(M,length(new.cases))] # get rid of unneeded cases (might be overly conservative)

        new.max <- max(new.cases,length(day.counts)) # maximal day in the data
        day.counts <- c(day.counts,rep(0,new.max-length(day.counts))) # 
        new.counts <- as.vector(table(new.cases))
        new.cts.tab <- rep(0,new.max)
        new.days <- unique(new.cases)
        new.cts.tab[sort(new.days)] <- new.counts
        
        day.counts <- day.counts+new.cts.tab  #allows to keep going if not close to the end

        min.ind <- min(which(new.cts.tab != 0))  
    }

    ind.max.dc <- max(which(cumsum(day.counts)<=(M+bdd)))
    sum.dc <- sum(day.counts[1:ind.max.dc])
    if((sum.dc <= (M+bdd)) & (sum.dc >= (M-bdd)))    day.counts <- day.counts[1:ind.max.dc]

    if((sum.dc <= (M+bdd)) & (sum.dc >= (M-bdd))){
        l.dc <- length(day.counts)
    
        if(m==1){
            Ns <- matrix(0,nr=1,nc=l.dc)
            Ns[m,] <- day.counts
        }
        if(m > 1){
            l.Ns <- length(Ns[m-1,])
            if(l.Ns > l.dc)  Ns <- rbind(Ns,c(day.counts,rep(0,l.Ns-l.dc)))
            if(l.Ns < l.dc)  Ns <- rbind(cbind(Ns,matrix(0,nc=(l.dc-l.Ns),nr=(m-1))),day.counts)
            if(l.Ns == l.dc) Ns <- rbind(Ns,day.counts)
        }
        
        m <- m+1
    }
    ind <- ind+1
 #   print(cbind(ind,m))
}
row.names(Ns) <- NULL

return(Ns)

}

####################################################################
## SIMULATE EPIDEMICS OF THE SAME LENGTH, SIZE IS NOT CONSTRAINED ##
####################################################################
# 11. MULTINOMIAL DISTRIBUTION #
getNsMNL <- function(R0,ps,num.samples=1000,max.days,num.start){
    # function to generate data for given parameter values
    
set.seed(21)

num.iter <- max.days# number of generations

R0 <- R0 
p <- ps 
k <- length(p)

# space allocation  
Ns <- matrix(1,nr=num.samples,nc=max.days)

for(m in 1:num.samples){
    day.cases <- rep(1,num.start)
    new.cases <- day.cases
    for(i in 1:num.iter){
        num.new <- rpois(length(new.cases),R0) # vector of the number of new cases generated by each case

        new.cases.start <- rep(new.cases,num.new)
        n <- sum(num.new) #total number of new cases on this iteration
        new.cases <- new.cases.start + c(1:k)%*%(rmultinom(n,1,p)) #cases that are new for this iteration

        day.cases <- c(day.cases,new.cases) #all cases
    }

    unique.days <- sort(unique(day.cases))
    day.counts <- as.vector(table(day.cases))

    N1 <- rep(0,max(unique.days))
    N1[unique.days] <- day.counts

    N <- N1[1:max.days]
    Ns[m,] <- N

    rm(list=c("num.new","day.cases","new.cases","unique.days","day.counts"))
}

return(Ns)

}

# 12. GAMMA DISTRIBUTION #
getNsGamL <- function(R0,shape,rate,num.samples=1000,max.days,num.start){
	# function to generate data for given parameter values
	# Length of the pidemic is constrained, size is not	
    
set.seed(21)

num.iter <- max.days# number of generations
 
# space allocation  
Ns <- matrix(1,nr=num.samples,nc=max.days)

for(m in 1:num.samples){
    day.cases <- rep(1,num.start)
    new.cases <- day.cases
    for(i in 1:num.iter){
        num.new <- rpois(length(new.cases),R0) # vector of the number of new cases generated by each case

        new.cases.start <- rep(new.cases,num.new)
        n <- sum(num.new) #total number of new cases on this iteration
        new.cases <- new.cases.start + ceiling(rgamma(n,shape=shape,rate=rate)) #cases that are new for this iteration

        day.cases <- c(day.cases,new.cases) #all cases
    }

    unique.days <- sort(unique(day.cases))
    day.counts <- as.vector(table(day.cases))

    N1 <- rep(0,max(unique.days))
    N1[unique.days] <- day.counts

    N <- N1[1:max.days]
    Ns[m,] <- N

    rm(list=c("num.new","day.cases","new.cases","unique.days","day.counts"))
}

return(Ns)

}


########################
# OTHER MISC FUNCTIONS #
########################
# 13. Getting summary stats for the serial interval #
getSummaryMN <- function(pars){
    # pars is in the form of (R0,ps)
    
    k <- length(pars)-1
    ps <- pars[2:(1+k)]
    mean.si <- sum(c(1:k)*ps)
    var.si <- sum((c(1:k)^2)*ps)-mean.si^2
    sd.si <- sqrt(var.si)   
    
    return(c(mean.si,sd.si))
}

###########################
# 14.  Function to get coverage probabilities for estimates of R0 and SI 
#      with WP method
#############################
covProbsWP <- function(parsEsts,k,epi.size,num.start,numRep=1000,distn="MN",est.si=T){
    
    r0 <- parsEsts[1]
    
    # get epis with the same number of cases-
    if(distn=="MN"){
        ps <- parsEsts[2:(k+1)]
        Ns <- getNsMN(r0,ps,num.samples=numRep,epi.size=epi.size,num.start=num.start)
    }
    if(distn=="gamma"){
        shape <- parsEsts[2]
        rate <- parsEsts[3]
        Ns <- getNsGam(r0,shape,rate,num.samples=numRep,epi.size=epi.size,num.start=num.start)
    }
        
    # perform estimation on the simulated epidemics
    if(est.si){
    tmp.ests <- matrix(0,nr=numRep,nc=length(parsEsts)+2) # thetas, len of epi, converge
    for(i in 1:numRep){
        Ns.tmp <- Ns[i,]
        Ns.tmp <- Ns.tmp[1:max(which(Ns.tmp!=0))]
        
        if(distn=="MN"){
            tmp <- constrOptim(rep(1.2,k),likN.MN,N=Ns.tmp,ui=diag(1,k),ci=rep(0,k),method="Nelder-Mead") # gets theta estimates
			pars.tmp <- c(sum(tmp$par),tmp$par/sum(tmp$par))
        }

        if(distn=="gamma"){
            tmp <- constrOptim(c(1.2,2,2),likN.gam,N=Ns.tmp,k=k,ui=diag(1,3),ci=rep(0,3),method="Nelder-Mead") # gets par estimates
			pars.tmp <- tmp$par
        }
		
		tmp.ests[i,] <- c(pars.tmp,length(Ns.tmp),tmp$convergence)
    } #end for loop

    covProbs <- rbind(apply(tmp.ests[,1:length(pars.tmp)],2,quantile,probs=c(0.025,0.5,0.975)),apply(tmp.ests[,1:length(pars.tmp)],2,min),apply(tmp.ests[,1:length(pars.tmp)],2,max))
	row.names(covProbs) <- c("2.5%","median","97.5%","min","max")

	} #end if(est.si)

    if(!est.si){
    tmp.ests <- matrix(0,nr=numRep,nc=2) # R0, len of epi
    for(i in 1:numRep){
        Ns.tmp <- Ns[i,]
        Ns.tmp <- Ns.tmp[1:max(which(Ns.tmp!=0))]
        
        	if(distn=="gamma"){
        		pars.tmp <- getR0.WP(N=Ns.tmp,pars=c(shape,rate),distn="gamma",k=20)			}
        	if(distn=="MN"){
        		pars.tmp <- getR0.WP(N=Ns.tmp,pars=ps,distn="MN",k=length(ps))        
        	}
			tmp.ests[i,] <- c(pars.tmp,length(Ns.tmp))
    } #end for loop
    covProbs <- c(quantile(tmp.ests[,1],probs=c(0.025,0.5,0.975)),min(tmp.ests[,1]),max(tmp.ests[,1]))
	names(covProbs) <- c("2.5%","median","97.5%","min","max")

	} #end if(!est.si)
    
    return(covProbs)
    
#    return(tmp.ests)
}

