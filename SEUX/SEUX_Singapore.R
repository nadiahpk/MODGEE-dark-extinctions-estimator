rm(list=ls())

library('seux')

source('get_CI_estimate.R')

# Compute overall fraction of extinctions including undetected extinctions.
# Arguments:
# * S_t = time series of detected, extant species (must start at zero)
# * E_t = time series of detected, extinct species (must start at zero)
# Return value:
# * p_t = time series of overall fraction of extinctions
compute_p_t <- function(S_t,E_t)
  {
  stopifnot(length(S_t)==length(E_t))
  # Method assumes that initially zero species have been detected
  stopifnot(S_t[1]==0 & E_t[1]==0)
  TT = length(S_t)-1
  mu_t_est = diff(E_t)/S_t[-(TT+1)]
  # By assumption, the extinction rate in the initial year is zero
  mu_t_est[1] = 0
  p_t = c(0,1-cumprod(1-mu_t_est))
  return(p_t)
  }

# Compute undetected extinctions assuming that the extinction rates for
# undetected species are the same as those for detected species.
# Arguments:
# * S_t = time series of detected, extant species (must start at zero)
# * E_t = time series of detected, extinct species (must start at zero)
# Return value:
# * X_t_est = time series of estimated undetected, extinct species
compute_X_t <- function(S_t,E_t)
  {
  p_t = compute_p_t(S_t,E_t)
  TT = length(S_t)-1
  X_t_est = S_t[TT+1]*p_t/(1-p_t[TT+1])-E_t
  return(X_t_est)
  }

# Compute undetected extinctions allowing for different extinction rates of
# detected and undetected species.
# Arguments:
# * S_t = time series of detected, extant species (must start at zero)
# * E_t = time series of detected, extinct species (must start at zero)
# * mu_t_star = time series of extinction rates for undetected species
# Return value:
# * X_t_est = time series of estimated undetected, extinct species
# The length of mu_t_star should be one less than the length of S_t and E_t,
# because it is assumed there is no information about mu_t_star at the final
# time (t=T).
# Gives the same results as compute_X_t() in the following special case:
# mu_t_est = c(diff(E_t)/S_t[-(TT+1)])
# mu_t_est[1] = 0
# X_t_est = compute_X_t_star(S_t,E_t,mu_t_est)
compute_X_t_star <- function(S_t,E_t,mu_t_star)
  {
  stopifnot(length(S_t)==length(E_t) & length(mu_t_star)==length(S_t)-1)
  
  TT0 = length(S_t)-1
  Y0 = (S_t[-1]+E_t[-1]-S_t[-(TT+1)]-E_t[-(TT+1)])/(E_t[TT+1]+S_t[TT+1])

  iii = which(S_t+E_t==max(S_t+E_t))
  # This assumes that the number of undetected species reached zero after the
  # last detection. This is necessary, because it's impossible to estimate
  # the rate of detection if no species are being detected (and so we
  # effectively assume the rate is zero after this time).
  TT = iii[1]-1

  Y = Y0[1:TT]
  mu_t_star = c(mu_t_star,0)

  temp = numeric(TT)
  
  for ( j in 0:(TT-2) )
    {
    temp[j+1] = prod(1-mu_t_star[((j+1):(TT-1))+1])
    }
  # Special case for j = TT-1 (in the product above, the start index would be
  # greater than the end index for this case, which means the product should be
  # equal to 1)
  temp[(TT-1)+1] = 1

  nu_est = numeric(TT+1)
  for ( i in 0:(TT-1) )
    {
    nu_est[i+1] = Y[i+1]*prod(1-mu_t_star[(i:(TT-1))+1])/sum((Y*temp)[(i:(TT-1))+1])
    }

  temp2 = numeric(TT-1+1)
  for ( j in 0:(TT-1) )
    {
    if ( j==0 )
      temp2[j+1] = mu_t_star[j+1]
    else
      temp2[j+1] = mu_t_star[j+1]*prod((1-mu_t_star[1:TT]-nu_est[1:TT])[(0:(j-1))+1])
    }
  
  X_t_est = (E_t[TT+1]+S_t[TT+1])*c(0,cumsum(temp2))/(1-sum(temp2))
  if ( TT!=TT0 ) X_t_est = c(X_t_est,rep(tail(X_t_est,1),TT0-TT))
  return(X_t_est)
  }

all_taxonomic_groups = c('mammals','bees','amphibians','reptiles','phasmids','crustaceans','fishes')
n = length(all_taxonomic_groups)
all_S_Ts = numeric(n)
all_E_Ts = numeric(n)
all_X_Ts = numeric(n)
all_p_Ts = numeric(n)
all_p_T_los = numeric(n)
all_p_T_his = numeric(n)

all_results = list()

########
graphics.off()
par(mfrow=c(1,1),mar=c(4,11,1,1))

draw_old_X = F

for ( group_i in 1:n )
  {
  taxonomic_group = all_taxonomic_groups[group_i]
  cat('**** Fitting to ',taxonomic_group,'\n',sep='')
  
  if ( taxonomic_group == 'mammals' )
    {
    cutoff_year = 2016
    } else if ( taxonomic_group == 'bees' )
    {
    cutoff_year = 2000
    } else if ( taxonomic_group == 'amphibians' )
    {
    cutoff_year = 2017
    } else if ( taxonomic_group == 'reptiles' )
    {
    cutoff_year = 1988
    } else if ( taxonomic_group == 'phasmids' )
    {
    cutoff_year = 2017
    } else if ( taxonomic_group == 'crustaceans' )
    {
    cutoff_year = 2007
    } else if ( taxonomic_group == 'fishes' )
    {
    cutoff_year = 2017
    }
  
  # Load in the data
  filename = paste(taxonomic_group,'_data.csv',sep='')
  mydata = read.csv(filename,header=T)
  
  stopifnot(all(mydata$First.record <= mydata$Last.record))
  
  year_min = min(mydata$First.record)
  year_max = max(mydata$Last.record)
  
  years = year_min:year_max
  
  # we assume that species seen in cutoff_year or after are extant
  mydata$Last.record[mydata$Last.record>=cutoff_year] = year_max
  
  # Let the first year be year t=1 and the last year be year t=TT.
  # So then t=0 is the year before the first year, where we have no data but by assumption mu=0.
  TT = length(years)
  
  S_t = numeric(TT+1)
  E_t = numeric(TT+1)
  
  # Compute the cumulative number of detected extant and detected extinct species in each year
  for ( t in 1:TT )
    {
    year = years[t]
    S_t[t+1] = sum(mydata$First.record<=year & mydata$Last.record>=year)
    E_t[t+1] = sum(mydata$Last.record<year)
    }
  
  # Estimate the number of undetected extinctions
  X_t_est = compute_X_t(S_t,E_t)
  
  # Estimate number of undetected extant species
  NN = S_t[TT+1]+E_t[TT+1]+X_t_est[TT+1]
  U_t_est = NN-S_t-E_t-X_t_est
  U_t_est_rounded = round(U_t_est) 
  
  # Bootstrap to get CI on the number of undetected extinct species
  n_rep = 1000
  all_X = matrix(NA,n_rep,length(X_t_est))
  
  mu_t_est = diff(E_t)/S_t[-(TT+1)]
  mu_t_est[1] = 0
  
  for ( j in 1:n_rep ){
    # This randomisation procedure generates new mu values from the binomial
    # distribution, based on the observed mu and S, and estimated U.
    # A resampled set of mu values is discarded if one of the mu values is 1, since this causes 
    # the corresponding estimate of X to be undefined.
    success_flag = 0
    while(success_flag == 0){
      mu_resample = rbinom(TT,S_t[-(TT+1)],mu_t_est)/S_t[-(TT+1)]
      mu_resample[1] = 0
      mu_star_resample = rbinom(TT,U_t_est_rounded[-(TT+1)],mu_resample)/U_t_est_rounded[-(TT+1)]
      U_t_est_rounded_indzero = which(U_t_est_rounded[-(TT+1)]==0)
      mu_star_resample[U_t_est_rounded_indzero] = 0
      if(max(mu_star_resample)<1){
        XX = compute_X_t_star(S_t,E_t,mu_star_resample)
        success_flag = 1
      }
    }
    all_X[j,] = XX
  }
  
  X_lo = apply(all_X,2,function(x) quantile(x,0.025))
  X_hi = apply(all_X,2,function(x) quantile(x,0.975))

  ####
  do_CI_est = ( E_t[TT+1] > 0 )

  if ( do_CI_est )
    {
    model_inputs <- get_model_inputs( mydata$First.record, mydata$Last.record )
    temp <- get_CI_estimate( model_inputs$S, model_inputs$E, nreps = 100000 )
    CIs_estimates <- temp[[1]]

    # organise X versus year data so that we can work with it in the old way
    years2a = c(head(rep(model_inputs$year,each=2)+rep(c(0,1),length(model_inputs$year)),-1),tail(years,1))
    # find duplicate years
    iii = which(diff(years2a)==0)
    if ( length(iii) ) years2 = years2a[-iii] else years2 = years2a
    kkka = c(tail(rep(1:length(model_inputs$year),each=2),-1),length(model_inputs$year))
    if ( length(iii) ) kkk = kkka[-iii] else kkk = kkka
    X_mean2 = CIs_estimates$X_mean[kkk]
    X_t_est2 = X_mean2
    X_lo2 = CIs_estimates$X_lo[kkk]
    X_hi2 = CIs_estimates$X_hi[kkk]
    XM2 = temp[[2]][,kkk]
    
    E_t2 = model_inputs$E[kkk]
    S_t2 = model_inputs$S[kkk]
    }

  # Plot the time series
  par(mar=c(5,5,2,1))
  yyy = S_t+E_t
  if ( do_CI_est ) yyy = c(yyy,X_hi2)
  ylim = c(0,max(yyy))
  plot(years,(S_t+E_t)[-1],type='l',ylim=ylim,
       xlab='Year',ylab='Number of species',lwd=2,cex.axis=2,cex.lab=2)
  lines(years,S_t[-1],col='green',lwd=2)
  lines(years,E_t[-1],col='red',lwd=2)

  if ( draw_old_X )
    {
    # Plot the time series of the number of undetected extinct species
    lines(years,X_t_est[-1],col='blue',lwd=2,lty=2)
    # Plot the 95% CI on undetected extinctions
    polygon(c(years,rev(years),years[1]),c(X_lo[-1],rev(X_hi[-1]),X_lo[2]),col=rgb(0,0,1,0.2),border=NA)
    }
    
  if ( do_CI_est )
    {
    lines(years2,X_mean2,col='magenta',lwd=2,lty=2)
    polygon(c(years2,rev(years2),years2[1]),
      c(X_lo2,rev(X_hi2),X_lo2[1]),col=rgb(1,0,1,0.2),border=NA)
    }
  
  ####
  
  if ( taxonomic_group == 'mammals' )
    legend(1895,20,
      legend=c(expression(paste('Total known (',italic(S[t]+E[t]),')',sep='')),
               expression(paste('Known extant (',italic(S[t]),')',sep='')),
               expression(paste('Known extinct (',italic(E[t]),')',sep='')),
               expression(paste('Inferred extinct (',italic(hat(X)[t]),') ',sep=''))
               ),
      col=c('black','green','red','magenta'),lty=c(1,1,1,2),lwd=2,
      cex=1)
  
  {
  # Summarise results
  cat('Number of detected extant species (S_T) = ',S_t[TT+1],'\n',sep='')
  cat('Number of detected extinct species (E_T) = ',E_t[TT+1],'\n',sep='')
  cat('Estimated overall proportion of extinctions (p_T) = ',
      round(tail(((E_t+X_t_est)/(S_t+E_t+X_t_est)),1),3),
      ' [',round(tail(((E_t+X_lo)/(S_t+E_t+X_lo)),1),3),',',
      round(tail(((E_t+X_hi)/(S_t+E_t+X_hi)),1),3),']\n',sep='')
  cat('Estimated number of undetected extinctions [assuming U_T=0] (X_T) = ',round(tail(X_t_est,1),1),
      ' [',round(tail(X_lo,1),1),',',round(tail(X_hi,1),1),']\n',sep='')
  }
  
  ####
  if ( do_CI_est )
    {
    cat('== with new CI estimates:\n')
    cat('Estimated overall proportion of extinctions (p_T) = ',
        round(tail(((E_t2+X_t_est2)/(S_t2+E_t2+X_t_est2)),1),3),
        ' [',round(tail(((E_t2+X_lo2)/(S_t2+E_t2+X_lo2)),1),3),',',
        round(tail(((E_t2+X_hi2)/(S_t2+E_t2+X_hi2)),1),3),']\n',sep='')
    cat('Estimated number of undetected extinctions [assuming U_T=0] (X_T) = ',round(tail(X_t_est2,1),1),
        ' [',round(tail(X_lo2,1),1),',',round(tail(X_hi2,1),1),']\n',sep='')
    
    all_S_Ts[group_i] = tail(S_t2,1)
    all_E_Ts[group_i] = tail(E_t2,1)
    all_X_Ts[group_i] = tail(X_t_est2,1)
    all_p_Ts[group_i] = tail((E_t2+X_t_est2)/(S_t2+E_t2+X_t_est2),1)
    all_p_T_los[group_i] = tail((E_t2+X_lo2)/(S_t2+E_t2+X_lo2),1)
    all_p_T_his[group_i] = tail((E_t2+X_hi2)/(S_t2+E_t2+X_hi2),1)
    }
  
  if ( do_CI_est )
    all_results[[group_i]] = list(years=years2,S_t=S_t2, E_t=E_t2, X_t_est=X_t_est2, XM=XM2) else
    all_results[[group_i]] = list(years=years,S_t=S_t, E_t=E_t, X_t_est=X_t_est, XM=matrix(0,100000,length(S_t)))
  names(all_results)[group_i] = taxonomic_group
  }

