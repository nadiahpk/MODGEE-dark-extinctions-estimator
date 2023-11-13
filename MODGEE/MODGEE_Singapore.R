rm(list=ls())

library(optimx) # for optimx()
library(compiler) # for cmpfun()
library(Rcpp)

# change this depending on which taxonomic group is wanted
taxonomic_group = 'butterflies'

if ( taxonomic_group == 'birds' )
  {
  yr_interval = 10
  n_detect_first_t = 10
  do_k = F; k_scale = 1
  do_l = T
  } else if ( taxonomic_group == 'butterflies' )
  {
  yr_interval = 10
  n_detect_first_t = 10
  # if k_scale=1 for butterflies the model doesn't converge
  do_k = F; k_scale = 10
  do_l = T
  } else if ( taxonomic_group == 'plants' )
  {
  yr_interval = 10
  n_detect_first_t = 10
  do_k = F; k_scale = 1
  do_l = F
  stopifnot(!do_l)
  }

Sys.setenv(PKG_CXXFLAGS = paste('-O3',Sys.getenv('PKG_CXXFLAGS'),sep=' '))
Rcpp::sourceCpp("compute_L_i_no_extinction.cpp")
Rcpp::sourceCpp("compute_L_i.cpp")

compute_L_i_no_extinction <- compute_L_i_no_extinction_Cpp
compute_L_i <- compute_L_i_Cpp

# compute the log-likelihood of the data given detection parameters
compute_neg_LL_detect_param <- function(parms,l_M,l_t_hats)
  {
  last_parms <<- parms
  l_n =  dim(l_M)[1]
  l_TT = dim(l_M)[2]
  # the first n-1 parameters are the detection rates for species [2,n]
  # (the detection rate for the first species is set to 1.0, without loss of generality)
  # the next T-1 parameters are the detection rates in years [1, T-1]
  # (the detection rate in year T can't be estimated using this method)
  # but note the matrix passed in as l_M here will already have been reduced to l_TT = T-1 columns
  l_ds = c(1.0,exp(parms[1:(l_n-1)]))
  l_hs = exp(parms[(1:l_TT)+l_n-1])
  Ls = lapply(1:l_n,function(i)
    {
    # t_max will be the year before the final year in which this species was detected
    t_max = l_t_hats[i]
    # compute the overall detection parameters for this species in years [1, t_max]
    xs_i = l_ds[i]*l_hs[1:t_max]
    # this is the likelihood of the parameters given M_i and no extinction
    L_i = compute_L_i_no_extinction(l_M[i,1:t_max],xs_i)
    return(L_i)
    })
  
  LL = sum(log(unlist(Ls)))
  stopifnot(!is.nan(LL))

  return(-LL)
  }

# compute the log-likelihood of the data given extinction parameters (detection parameters fixed)
compute_neg_LL_extinct_param <- function(parms,l_M,l_H,l_ds,l_hs,do_k=T,do_l=T)
  {
  last_parms <<- parms
  l_n =  dim(l_M)[1]
  l_TT = dim(l_M)[2]

  if ( do_l) l_n_H = dim(l_H)[2]
    
  # the first parameter is the detection rate in the final timestep T
  # (the detection rates in timesteps 1 to T-1 are passed in as fixed arguments)
  l_h_T = exp(unlist(parms[1]))
  # the next T-2 parameters are the extinction rates in years [3,T]
  # (the extinction rate in the first year can't be estimated because there are no prior observations)
  # (the extinction rate in the second year can't be estimated either because we are censoring
  #  rows with only one observation, so there will be no {1,0,0,...0} rows)
  # (the extinction rate in the final year can't be estimated because it can't be untangled from the
  #  detection rate in that year)
  l_ms = c(0,0,unlist(exp(parms[1+(1:(l_TT-3))])),0)
  if ( do_k ) l_k=parms[1+(l_TT-3)+1]/k_scale else l_k=0
  if ( do_l ) l_ls=c(0,unlist(parms[1+(l_TT-3)+do_k+(1:(l_n_H-1))])) else l_ls=NULL

  zeros = rep(0,l_TT)
  ones = lapply(1:l_TT,function(i) {x=zeros;x[i]=1;return(x)})
  Ls = lapply(1:l_n,function(i)
    {
    #l_ms_prime = l_ds[i]^l_k*l_ms*exp(sum(l_H[i,]*l_ls))
    l_ms_prime = l_ms
    if ( do_k ) l_ms_prime = l_ds[i]^l_k*l_ms_prime
    if ( do_l ) l_ms_prime = l_ms_prime*exp(sum(l_H[i,]*l_ls))
    # this is the probability observing M_i given the parameters
    P_i_unconditional = compute_L_i(l_M[i,],l_ds[i]*c(l_hs,l_h_T),l_ms_prime)
    # this is the probability of observing M_i = 0 given the parameters
    P_i_M_zero = compute_L_i(zeros,l_ds[i]*c(l_hs,l_h_T),l_ms_prime)
    # this is the probability of observing M_i containing only one instance of 1 given the parameters
    P_i_M_one = unlist(lapply(ones,function(ones_i) compute_L_i(ones_i,l_ds[i]*c(l_hs,l_h_T),l_ms_prime)))
    # this is the conditional likelihood given that M_i has at least two instances of 1
    L_i_conditional = P_i_unconditional/(1-P_i_M_zero-sum(P_i_M_one))
    return(L_i_conditional)
    })
  LL = sum(log(unlist(Ls)))
  stopifnot(!is.nan(LL))

  return(-LL)
  }

# compute the log-likelihood of the data given the detection parameter for a species that was recorded
# only once (all other parameters fixed)
compute_neg_LL_detect_param_1_obs <- function(parm,l_M,l_H,l_ms,l_hs,l_k,l_ls,do_k=T,do_l=T)
  {
  last_parm <<- parm
  l_TT = length(l_M)

  l_d = exp(parm)

  zeros = rep(0,l_TT)
  ones = lapply(1:l_TT,function(i) {x=zeros;x[i]=1;return(x)})
  l_ms_prime = l_ms
  if ( do_k ) l_ms_prime = l_d^l_k*l_ms_prime
  if ( do_l ) l_ms_prime = l_ms_prime*exp(sum(l_H*l_ls))
  # this is the probability of observing M_i given the parameter
  P_i_unconditional = compute_L_i(l_M,l_d*l_hs,l_ms_prime)
  # this is the probability of observing M_i containing only one instance of 1 given the parameter
  P_i_M_one = unlist(lapply(ones,function(ones_i) compute_L_i(ones_i,l_d*l_hs,l_ms_prime)))
  # this is the conditional likelihood of the parameter given that M_i contains only one instance of 1
  LL = log(P_i_unconditional/sum(P_i_M_one))
  stopifnot(!is.nan(LL))

  return(-LL)
  }

# fit the detection and extinction parameters to an observation matrix using the two-step method
est_detect_extinct_parms_ML <- function(M,H,ds_start,hs_start,ms_start,k_start=NULL,ls_start=NULL)
  {
  n_spp = dim(M)[1]
  TT = dim(M)[2]

  do_k = !is.null(k_start)
  do_l = !is.null(ls_start)

  if ( do_l) n_H = dim(H)[2]
      
  cat('do_k = ',do_k,'\n',sep='')
  cat('do_l = ',do_l,'\n',sep='')
  
  # M_prime will be a binary matrix with 1 at entry {i,j} indicating that species i was
  # definitely extant at time j+1 (because it was detected at time j+1 or any time thereafter)
  M_prime = t(apply(M,1,function(temp)
    {
    temp2 = 0*temp
    iii = which(temp==1)
    if ( length(iii) ) if ( max(iii)>1 ) temp2[1:(max(iii)-1)] = 1
    return(temp2)
    }))
  
  n_obs = rowSums(M)
  # take only species with at least two observations
  # because otherwise we can't estimate the detection rate
  iii = which(n_obs>=2)
  not_iii = which(n_obs<2)
  n_fit = length(iii)
  # t_hats[i] gives the time just prior to the final detection of species i
  t_hats = apply(M_prime,1,function(M_prime_i) {x=tail(which(M_prime_i==1),1);ifelse(length(x),x,NA)})

  # the detection rate of the first species is fixed at 1.0
  # the detection rate in the final timestep can't be estimated in the first step
  parms = c(as.numeric(log(ds_start[iii][-1])),as.numeric(log(hs_start[-TT])))
  
  cat('Fitting detection parameters...\n')
  best_parms = try(optimx(parms,
                      compute_neg_LL_detect_param,l_M=M[iii,-TT],l_t_hats=t_hats[iii],
                      method="BFGS"))
  
  if ( 'try-error'%in%class(best_parms) | length(best_parms)==1 | any(is.na(best_parms[1:length(parms)])) )
    {
    cat('Could not fit detection parameters!\n')
    return(NULL)
    }
   
  ds_est = unlist(c(1.0,exp(best_parms[1:(n_fit-1)])))
  hs_est = unlist(c(exp(best_parms[n_fit-1+(1:(TT-1))]),NA))
  
  # the extinction rates in the first, second and final timesteps can't be estimated
  parms = c(log(hs_start[TT]),log(ms_start[-c(1,2,TT)]))
  if ( do_k ) parms=c(parms,k_start)
  if ( do_l ) parms=c(parms,ls_start[-1])

  cat('Fitting extinction parameters...\n')
  best_parms = try(optimx(parms,
                      compute_neg_LL_extinct_param,l_M=M[iii,],l_H=H[iii,],l_ds=ds_est,l_hs=hs_est[-TT],
                      do_k=do_k,do_l=do_l,
                      method="BFGS"))

  if ( 'try-error'%in%class(best_parms) | length(best_parms)==1 | any(is.na(best_parms[1:length(parms)])) )
    {
    cat('Could not fit extinction parameters!\n')
    return(NULL)
    }
  
  hs_est[TT] = unlist(exp(best_parms[1]))
  ms_est = c(0,0,exp(unlist(best_parms[1+(1:(TT-3))])),0)
  if ( do_k ) k_est = best_parms[1+(TT-3)+1] else k_est = NULL
  if ( do_l ) ls_est = c(0,unlist(best_parms[1+(TT-3)+do_k+(1:(n_H-1))])) else ls_est = NULL
  LL = -unlist(best_parms["value"])

  ds_est_final = rep(NA,n_spp)
  ds_est_final[iii] = ds_est

  small_d = 1e-5
  big_d = 100
  
  if ( length(not_iii) )
    {
    cat('Fitting detection parameters for species with only one detection...\n', sep='')
    for ( j in 1:length(not_iii) )
      {
      i = not_iii[j]
      n_guess=1000
      d_guesses = exp(seq(log(small_d),log(big_d),length.out=n_guess))
      neg_LL_guesses = rep(NA,n_guess)
      for ( k in 1:n_guess )
        {
        neg_LL_guesses[k] = compute_neg_LL_detect_param_1_obs(log(d_guesses[k]),
          l_M=M[i,],l_H=H[i,],l_ms=ms_est,l_hs=hs_est,l_k=k_est,l_l=ls_est,do_k=do_k,do_l=do_l)
        }
      plot(d_guesses,exp(-neg_LL_guesses),type='l',log='x')
      parm = log(d_guesses[which.min(neg_LL_guesses)])
      best_parm = try(optimx(parm,
                          compute_neg_LL_detect_param_1_obs,l_M=M[i,],l_H=H[i,],
                          l_ms=ms_est,l_hs=hs_est,l_k=k_est,l_l=ls_est,
                          do_k=do_k,do_l=do_l,
                          method="BFGS"))
    
      if ( 'try-error'%in%class(best_parm) | is.na(best_parm[[1]]) )
        {
        #cat('Could not fit detection parameter for species ',i,'!\n',sep='')
        d_est_final[i] = d_guesses[which.max(-neg_LL_guesses)]
        } else
        {
        ds_est_final[i] = exp(best_parm[[1]])
        }
      }
    
    # don't allow the estimated d to be too big, because otherwise it causes numerical errors later
    big_d = 100
    iii_too_big = which(ds_est_final[not_iii]>big_d)
    ds_est_final[not_iii][iii_too_big] = big_d
    }
  
  r_val = list(M=M,H=H,ds_est=ds_est_final,hs_est=hs_est,ms_est=ms_est,k_est=k_est,ls_est=ls_est)
  
  r_val$LL = LL
  return(r_val);
  }

est_detect_parms <- function(M)
  {
  # M is spp_detected
  # M_prime is spp_extant_hat

  M_prime = t(apply(M,1,function(temp)
    {
    temp2 = 0*temp
    iii = which(temp==1)
    if ( length(iii) ) if ( max(iii)>1 ) temp2[1:(max(iii)-1)] = 1
    return(temp2)
    }))
   
  TT = dim(M)[2]
  stopifnot(dim(M_prime)[2]==TT)

  spp_detectability_hat = rowSums(M & M_prime)/rowSums(M_prime)

  delta_hat_new_t = numeric(TT)
  mean_detectability_hat = mean(spp_detectability_hat,na.rm=T)

  for ( t in 1:TT )
    {
    k = sum(M_prime[,t]) # note the k ends up cancelling, but I put it here for clarity
    delta_hat_t = sum(M[,t] & M_prime[,t])/k
    mu_t = sum(M[M_prime[,t]==1])/k
    delta_hat_new_t[t] = mean_detectability_hat*delta_hat_t/mu_t
    }

  d_est0 = spp_detectability_hat/mean(delta_hat_new_t/0.5,na.rm=T)
  h_est0 = delta_hat_new_t/0.5

  d_est = d_est0/d_est0[1]
  h_est = h_est0*d_est0[1]

  iii = c(TT,which(h_est==0))
  # detection effort in last year can't be estimated using the above method
  h_est[iii] = mean(h_est[-iii])
  
  # for species that were only seen once, d_est will be zero, but we set it
  # instead to NA, because it can't really be estimated
  iii = which(d_est==0)
  d_est[iii] = NA
  
  return(list(d_est=d_est,h_est=h_est))
  }

# create fake data given parameter values
rdetect <- function(ds, hs, ms)
  {
  n = length(ds)
  TT = length(hs)
  stopifnot(length(hs)==length(ms))
  mus = 1-exp(-ms)
  
  iii = 1:n

  M = matrix(NA,n,TT)
  M_extant = matrix(NA,n,TT)

  # the last value here will be the probability of the species still being extant after the final timestep
  p_extinction = c(mus,1)*c(1,cumprod(1-mus))
  
  n_attempt = 0

  while (length(iii))
    {
    n_attempt = n_attempt + 1
    M0 = matrix(runif(length(iii)*TT)<1-exp(-ds[iii]%o%hs),length(iii),TT)

    t_extinct = sample(1:(TT+1),length(iii),prob=p_extinction,replace=T)

    M1 = 1*do.call('rbind',lapply(1:length(iii), function(i)
      {
      return(c(rep(1,t_extinct[i]-1),rep(0,TT-t_extinct[i]+1)))
      }))

    M_extant[iii,] = M1
    M[iii,] = M0*M1

    iii = which(rowSums(M)<=1)
    }

  return(list(M_extant=M_extant,M=M))
  }

rdetect_multiple_subgroups <- function(ds, hs, ms_prime, ms_indices)
  {
  stopifnot(length(ds)==length(ms_indices))
  
  k = nrow(ms_prime)
  M = matrix(NA,length(ds),length(hs))
  M_extant = matrix(NA,length(ds),length(hs))
  for ( i in 1:k )
    {
    iii = which(ms_indices==i)
    if ( length(iii) ) 
      {
      temp_i = rdetect(ds[iii],hs,ms_prime[i,])
      M[iii,] = temp_i[[2]]
      M_extant[iii,] = temp_i[[1]]
      }
    }
  return(list(M_extant=M_extant,M=M))
  }

compute_P_extinct_i <- function(M_i,d,hs,ms)
  {
  E_t_ext_num = 0
  E_t2_ext_num = 0
  E_t_ext_denom = 0
  TT = length(M_i)
  ps = rep(0,TT)
  deltas = 1-exp(-d*hs)
  mus = 1-exp(-ms)
  qs = M_i*deltas+(1-M_i)*(1-deltas)
  log_p_M_i = log(compute_L_i(M_i,d*hs,ms))
  for ( t in 2:TT )
    {
    log_p_M_i_given_ext_t = sum(log(qs[1:(t-1)]))+sum(log((1-M_i)[t:TT]))
    log_p_ext_t = sum(log((1-mus)[1:(t-1)]))+log(mus[t])
    log_p_ext_t_given_M_i = log_p_M_i_given_ext_t+log_p_ext_t-log_p_M_i
    p_ext_t_given_M_i = exp(log_p_ext_t_given_M_i)
    ps[t] = ps[t-1] + p_ext_t_given_M_i
    E_t_ext_num = E_t_ext_num + t*p_ext_t_given_M_i
    E_t2_ext_num = E_t2_ext_num + t^2*p_ext_t_given_M_i
    E_t_ext_denom = E_t_ext_denom + p_ext_t_given_M_i
    }
  E_t_ext = E_t_ext_num/E_t_ext_denom
  E_t2_ext = E_t2_ext_num/E_t_ext_denom
  var_t_ext = E_t2_ext-E_t_ext^2
  
  return(list(ps,E_t_ext,var_t_ext))
  }

compute_summary_stats_subgroups <- function(l_fit)
  {
  stopifnot(do_l)
  
  l_M=l_fit$M
  l_H=l_fit$H
  l_ds=l_fit$ds
  l_hs=l_fit$hs
  l_ms=l_fit$ms
  l_ls=l_fit$ls

  l_n = nrow(l_M)
  TT = ncol(l_M)
  l_n_H = ncol(l_H)
  
  stopifnot(nrow(l_H)==l_n & length(l_ds)==l_n & length(l_hs)==TT & length(l_ms)==TT & length(l_ls)==l_n_H)

  p_exts = matrix(NA,l_n,TT)
  E_t_exts = rep(NA,l_n)
  var_t_exts = rep(NA,l_n)
  for ( j in 1:l_n )
    {
    ms_prime = l_ms*exp(sum(l_ls*l_H[j,]))
    temp = compute_P_extinct_i(l_M[j,],l_ds[j],l_hs,ms_prime)
    if ( !is.na(tail(temp[[1]],1)) )
      {
      p_exts[j,] = temp[[1]]
      E_t_exts[j] = temp[[2]]
      var_t_exts[j] = temp[[3]]
      }
    }

  H_list = list()
  for ( j in 1:l_n_H ) H_list[[j]] = c(0,1)
  unique_H = expand.grid(H_list)

  S_hat = t(apply(unique_H,1,function(unique_H_j) colSums(do.call('rbind',lapply(1:nrow(l_H),function(i)
    {
    if(all(unique_H_j==l_H[i,])) 1-p_exts[i,] else rep(NA,TT)
    }
    )),na.rm=T)))

  E_hat = t(apply(unique_H,1,function(unique_H_j) colSums(do.call('rbind',lapply(1:nrow(l_H),function(i)
    {
    if(all(unique_H_j==l_H[i,])) p_exts[i,] else rep(NA,TT)
    }
    )),na.rm=T)))
  
  ms_prime = t((l_ms%o%exp(as.matrix(unique_H)%*%l_ls))[,,1])
  Mus = t(apply(ms_prime,1,function(x) 1-cumprod(exp(-x))))

  E_plus_X_hat = S_hat[,TT]*Mus/(1-Mus[,TT])
  X_hat = E_plus_X_hat - E_hat

  # formula can be obtained by manipulating Eq. (10) in Chisholm et al. (2016)
  Mus_overall = colSums(E_plus_X_hat)/sum(S_hat[,TT]+E_plus_X_hat[,TT])

  return(list(Mus_overall=Mus_overall,p_exts=p_exts[,TT],E_t_exts=E_t_exts,var_t_exts=var_t_exts,
    S_T_hat=colSums(S_hat)[TT],E_T_hat=colSums(E_hat)[TT],p_T_hat=Mus_overall[TT],X_T_hat=colSums(X_hat)[TT],
    X_hat_by_subgroup=X_hat))
  }

compute_summary_stats <- function(l_fit)
  {
  if ( do_l ) return(compute_summary_stats_subgroups(l_fit))
  
  l_M=l_fit$M
  l_ds=l_fit$ds
  l_hs=l_fit$hs
  l_ms=l_fit$ms

  Mus_overall=1-cumprod(exp(-l_ms))

  n = dim(l_M)[1]
  p_exts = rep(NA,n)
  E_t_exts = rep(NA,n)
  var_t_exts = rep(NA,n)
  
  for ( i in 1:n )
    {
    if ( is.na(l_ds[i]) ) next
    temp = compute_P_extinct_i(l_M[i,],l_ds[i],l_hs,l_ms)
    p_exts[i] = tail(temp[[1]],1)
    E_t_exts[i] = temp[[2]]
    var_t_exts[i] = temp[[3]]
    }

  S_T_hat = sum(1-p_exts,na.rm=T)
  E_T_hat = sum(p_exts,na.rm=T)
  p_T_hat = Mus_overall[TT]
  X_T_hat = S_T_hat*p_T_hat/(1-p_T_hat)-E_T_hat
  
  return(list(Mus_overall=Mus_overall,p_exts=p_exts,E_t_exts=E_t_exts,var_t_exts=var_t_exts,
    S_T_hat=S_T_hat,E_T_hat=E_T_hat,p_T_hat=p_T_hat,X_T_hat=X_T_hat))
  }

# Load in the data
filename = paste(taxonomic_group,'_data.csv',sep='')
mydata = read.table(filename,sep=',',quote="",check.names=F,header=T)

yrs = as.numeric(names(mydata))
yr0 = min(yrs)
yr1 = max(yrs)
temp = floor(yrs/yr_interval)
ts = temp-min(temp)+1

####

spp_detected = t(apply(mydata,1,function(x) 1*by(x,ts,any)))

jjj = rev(order(rowSums(spp_detected)))
t0 = min(which(colSums(spp_detected)>=n_detect_first_t))
t1 = max(ts)
ttt = t0:t1
M = spp_detected[jjj,ttt]

if ( any(rowSums(M)==0) )
  {
  cat('Error: at least row of detection matrix is zero after trimming.')
  stopifnot(F)
  }

if (do_l)
  {
  if ( taxonomic_group=='birds' )
    {
    subgroup_data = read.table('Habitat_matrix_birds.csv',sep=',',quote="\"",header=T)
    
    h01 = subgroup_data$Forest.interior..Primary.forest.
    h02 = subgroup_data$Edge.forest..Secondary.forest.
    h03 = subgroup_data$Coastal.Mangrove
    h04 = subgroup_data$Grassland
    h05 = subgroup_data$Parkland
    
    h01[is.na(h01)] = 0
    h02[is.na(h02)] = 0
    h03[is.na(h03)] = 0
    h04[is.na(h04)] = 0
    h05[is.na(h05)] = 0
    
    h1 = ( h01 | h03 ) & !h02 & !h04 & !h05
    h2 = h02 & !h04 & !h05
    h3 = h04 | h05
    
    iii = match(rownames(mydata),subgroup_data$Common.name)
    
    spp_subgroup = cbind(h1,h2,h3)[iii,]
    
    H = spp_subgroup[jjj,]
    } else if ( taxonomic_group=='butterflies' )
    {
    family_data = read.table('butt_families.csv',sep=',',quote="\"",header=F)
    
    h1 = (family_data[,2]=='Lycaenidae')
    h2 = (family_data[,2]!='Lycaenidae')

    iii = match(rownames(mydata),family_data[,1])
    spp_subgroup = cbind(h1,h2)[iii,]
    H = spp_subgroup[jjj,]
    } else stopifnot(F)
  } else
  {
  H = NULL
  }

TT = dim(M)[2]

temp = est_detect_parms(M)
ds_start = rep(NA,nrow(M))
ds_start = temp$d_est
hs_start = temp$h_est
ms_start = c(0,0,runif(TT-2,0,1))

k_start = if (do_k) 0 else NULL
ls_start = if (do_l) rep(0,ncol(H)) else NULL

cat('Doing ML fit to ',taxonomic_group,' data...\n',sep='')
time0 = proc.time()
myfit = est_detect_extinct_parms_ML(M,H,ds_start,hs_start,ms_start,k_start=k_start,ls_start=ls_start)
time1 = proc.time()
print(time1-time0)

ds = myfit$ds_est
hs = myfit$hs_est
ms = myfit$ms_est
if ( do_k ) k = myfit$k_est/k_scale
if ( do_l ) ls = myfit$ls_est
LL = unlist(myfit$LL)
n_fit = sum(!is.na(ds))
n_parm = n_fit-1+TT+TT-3+do_k+do_l*length(ls)
AIC = 2*n_parm-2*LL

# don't do any subsequent analyses if do_K=TRUE, because this model is just used to see if its AIC
# is better than that of the do_K=FALSE model
stopifnot(!do_k)

####

fit_stats = compute_summary_stats(myfit)

yr_offset = min(floor(yrs/yr_interval)+(t0-1))
ts_yr = ((1:TT)-1+yr_offset)*yr_interval

####

graphics.off()
par(mar=c(4,4,1,1))

if ( !do_l )
  {
  cat('Cumulative extinction rates: ')
  print(round(fit_stats$Mus_overall,3))
  cat('\n')
  plot(ts_yr,fit_stats$Mus_overall,type='o',pch=19,ylim=range(c(0,fit_stats$Mus_overall)))
  lines(ts_yr,fit_stats$Mus_overall,type='o',pch=1,lty=2)
  }

if ( do_l )
  {
  compute_Mus <- function(ms,ls,l_H)
    {
    ms_prime = ms*exp(sum(ls*l_H))
    Mus = 1-cumprod(exp(-ms_prime))
    }

  if ( taxonomic_group == 'birds' )
    {
    Mus_primary = compute_Mus(ms,ls,c(1,0,0))
    Mus_secondary = compute_Mus(ms,ls,c(0,1,0))
    Mus_grassland_parkland = compute_Mus(ms,ls,c(0,0,1))
    par(mar=c(5,5,1,1))
    plot(ts_yr,Mus_primary,type='o',pch=19,ylim=range(c(0,1)),col='forest green',
      xlab='Decade',ylab='Estimated extinction rate',cex.axis=2,cex.lab=2)
    lines(ts_yr,Mus_secondary,type='o',pch=19,lty=1,col='brown')
    lines(ts_yr,Mus_grassland_parkland,type='o',pch=19,lty=1,col='pink')
    legend(1865,1.05,legend=c('Primary + old secondary','Young secondary','Grassland + parkland'),lty=1,pch=19,
      col=c('forest green','brown','pink'),cex=1.2)
    } else if ( taxonomic_group == 'butterflies' )
    {
    Mus_lycaenid = compute_Mus(ms,ls,c(1,0))
    Mus_other = compute_Mus(ms,ls,c(0,1))
    par(mar=c(5,5,1,1))
    plot(ts_yr,Mus_other,type='o',pch=19,ylim=range(c(0,0.7)),col='orange',
      xlab='Decade',ylab='Estimated extinction rate',cex.axis=2,cex.lab=2)
    lines(ts_yr,Mus_lycaenid,type='o',pch=19,lty=1,col='blue')
    legend(1850,0.7,legend=c('Lycaenids','Other'),lty=1,pch=19,
      col=c('blue','orange'),cex=1.2)
    }
  }

####
# calculate Pr(extinct) and E(t_extinct | extinct) for each species

####
# histogram of Pr(extinct)
graphics.off()
par(mar=c(4,5,2,2))
hist(fit_stats$p_exts,nclass=20,xlab='Pr(extinct)',ylab='Number of species',cex.axis=1.5,cex.lab=1.5,main='')

####
# graph of E(t_extinct | extinct) vs. Pr(extinct)
E_t_exts_decade = (fit_stats$E_t_exts-1+yr_offset)*yr_interval
plot(fit_stats$p_exts,E_t_exts_decade,pch=19,col=rgb(0.2,0.2,0.2),xlab='Pr(extinct)',ylab='Estimated date of extinction',
  cex.axis=1.5,cex.lab=1.5,log='x')

if ( taxonomic_group == 'plants' )
  {
  focal_spp = c('Homalomena griffithii','Micropera fuscolutea','Erycibe maingayi','Syzygium avene')
  lll = match(tolower(focal_spp),rownames(M))
  points(fit_stats$p_exts[lll],E_t_exts_decade[lll],pch=19,col='red')
  }

####
# histogram of conditional extinction dates
temp = hist(E_t_exts_decade,xlab='E(t_extinct | extinct)',ylab='Number of species',cex.axis=1.5,cex.lab=1.5,main='')

####
# check that extinction time is after the last detection date
#last_ts = apply(M,1,function(x) tail(which(x==1),1))
#last_ts_decade = last_ts*10+1860
#plot(last_ts_decade+10,E_t_exts_decade)
#abline(0,1,lty=2)

####
do_boot_tests = T
do_sim_tests = F

####
# do tests on bootstrapped data
if ( do_boot_tests )
  {
  n_rep = 200
  
  Mus_overall_boot = numeric(n_rep)
  
  all_boot_fits = list()
  all_boot_stats = list()
  
  time00 = proc.time()
  
  for ( i in 1:n_rep )
    {
    n = dim(M)[1]
    # sorting here ensures that the matrix rows are arranged in order
    # from most-detected to least-detected species again
    kkk = sort(sample(1:n,n,replace=T))
    M_boot = M[kkk,]
    H_boot = H[kkk,]
    temp = est_detect_parms(M_boot)
    ds_start = temp$d_est
    hs_start = temp$h_est
    ms_start = c(0,0,runif(TT-3,0,1),0)

    k_start = NULL
    ls_start = if(do_l) rep(0,ncol(H_boot)) else NULL
    
    cat('Doing ML fit to bootstrapped data (rep ',i,'/',n_rep,')...\n')
    time0 = proc.time()
    myfit_boot = est_detect_extinct_parms_ML(M_boot,H_boot,ds_start,hs_start,ms_start,k_start,ls_start)
    if ( is.null(myfit_boot) ) next
    time1 = proc.time()
    print(time1-time0)
    
    all_boot_fits[[i]] = myfit_boot
    if ( class(myfit_boot)=="NULL" ) all_boot_stats[[i]] = NULL else
      all_boot_stats[[i]] = compute_summary_stats(myfit_boot)
    }
  
  time01 = proc.time()
  print(time01-time00)

  all_Mus_overall_boot = do.call('rbind',lapply(all_boot_stats, function(x)
      {
      if ( class(x)=='NULL') return(rep(NA,TT)) else return(x$Mus_overall)
      }))

  Mus_boot_mean = apply(all_Mus_overall_boot,2,mean,na.rm=T)
  Mus_boot_lo = apply(all_Mus_overall_boot,2,quantile,0.025,na.rm=T)
  Mus_boot_hi = apply(all_Mus_overall_boot,2,quantile,0.975,na.rm=T)
  
  par(mar=c(5,5,1,1))
  plot(ts_yr,fit_stats$Mus_overall,type='o',pch=19,ylim=range(c(0,fit_stats$Mus_overall,Mus_boot_lo,Mus_boot_hi)),
    cex.lab=2,cex.axis=2,xlab='Decade',ylab='Estimated extinction rate')
  lines(ts_yr,Mus_boot_mean,type='o',pch=1,col='blue')
  polygon(c(ts_yr,rev(ts_yr),ts_yr[1]),c(Mus_boot_lo,rev(Mus_boot_hi),Mus_boot_lo[1]),col=rgb(0,0,1,0.3),border=NA)
  
  legend(1,0.8,legend=c('Fits to original data','Fits to boot-strapped data'),lty=1,pch=c(19,1),col=c('black','blue'),
    cex=1.5)
  }

####
# do tests on simulated data
if ( do_sim_tests )
  {
  n_rep = 200
  
  Mus_overall_sim = numeric(n_rep)
  
  all_sim_fits = list()
  all_sim_stats = list()
  
  time00 = proc.time()
  
  for ( i in 1:n_rep )
    {
    if ( do_l )
      {
      stopifnot(taxonomic_group=='birds')
      unique_H = expand.grid(h1=c(0,1),h2=c(0,1),h3=c(0,1))
      ms_prime = t((ms%o%exp(as.matrix(unique_H)%*%ls))[,,1])
      ms_indices = (H[,1]+H[,2]*2+H[,3]*4+1) # a bit of a hack
      temp = rdetect_multiple_subgroups(ds, hs, ms_prime, ms_indices)
      iii = which(!is.na(ds))
      H_sim = H[iii,]
      M_sim = temp[[2]][iii,]
      } else
      {
      temp = rdetect(ds, hs, ms)
      H_sim = NULL
      iii = which(!is.na(ds))
      M_sim = temp[[2]][iii,]
      }

    temp = est_detect_parms(M_sim)
    ds_start = temp$d_est
    hs_start = temp$h_est
    ms_start = c(0,0,runif(TT-3,0,1),0)

    k_start = NULL
    ls_start = if (do_l) rep(0,ncol(H)) else NULL
            
    cat('Doing ML fit to simulated data (rep ',i,'/',n_rep,')...\n')
    time0 = proc.time()
    myfit_sim = est_detect_extinct_parms_ML(M_sim,H_sim,ds_start,hs_start,ms_start,k_start,ls_start)
    if ( is.null(myfit_sim) ) next
    time1 = proc.time()
    print(time1-time0)
    
    all_sim_fits[[i]] = myfit_sim
    if ( class(myfit_sim)=="NULL" ) all_sim_stats[[i]] = NULL else
      all_sim_stats[[i]] = compute_summary_stats(myfit_sim)
    }
  
  time01 = proc.time()
  print(time01-time00)

  all_Mus_overall_sim = do.call('rbind',lapply(all_sim_stats, function(x)
      {
      if ( class(x)=='NULL') return(rep(NA,TT)) else return(x$Mus_overall)
      }))

  Mus_sim_mean = apply(all_Mus_overall_sim,2,mean,na.rm=T)
  Mus_sim_lo = apply(all_Mus_overall_sim,2,quantile,0.025,na.rm=T)
  Mus_sim_hi = apply(all_Mus_overall_sim,2,quantile,0.975,na.rm=T)
  
  par(mar=c(5,5,1,1))
  plot(ts_yr,fit_stats$Mus_overall,type='o',pch=19,ylim=range(c(0,fit_stats$Mus_overall,Mus_sim_lo,Mus_sim_hi)),
    cex.lab=2,cex.axis=2,xlab='Decade',ylab='Estimated extinction rate')
  lines(ts_yr,Mus_sim_mean,type='o',pch=1,col='red')
  polygon(c(ts_yr,rev(ts_yr),ts_yr[1]),c(Mus_sim_lo,rev(Mus_sim_hi),Mus_sim_lo[1]),col=rgb(1,0,0,0.3),border=NA)
  
  legend(1,0.8,legend=c('Fits to original data','Fits to simulated data'),lty=1,pch=c(19,1),col=c('black','red'),
    cex=1.5)
  }

####
{
cat('S_T_hat = ',fit_stats$S_T_hat,'\n',sep='')
cat('E_T_hat = ',fit_stats$E_T_hat,'\n',sep='')
cat('X_T_hat = ',round(fit_stats$X_T_hat,1),'\n',sep='')
cat('p_T_hat = ',round(fit_stats$p_T_hat,3),'\n',sep='')

if ( do_boot_tests )
  {
  cat('Boot-strapping:\n')
  S_T_hats = unlist(lapply(all_boot_stats, function(x) x$S_T_hat))
  E_T_hats = unlist(lapply(all_boot_stats, function(x) x$E_T_hat))
  X_T_hats = unlist(lapply(all_boot_stats, function(x) x$X_T_hat))
  p_T_hats = unlist(lapply(all_boot_stats, function(x) x$p_T_hat))
  cat('    S_T_hat = ',round(mean(S_T_hats),1),' [',round(quantile(S_T_hats,0.025),1),', ',round(quantile(S_T_hats,0.975),1),']\n',sep='')
  cat('    E_T_hat = ',round(mean(E_T_hats),1),' [',round(quantile(E_T_hats,0.025),1),', ',round(quantile(E_T_hats,0.975),1),']\n',sep='')
  cat('    X_T_hat = ',round(mean(X_T_hats),1),' [',round(quantile(X_T_hats,0.025),1),', ',round(quantile(X_T_hats,0.975),1),']\n',sep='')
  cat('    p_T_hat = ',round(mean(p_T_hats),3),' [',round(quantile(p_T_hats,0.025),3),', ',round(quantile(p_T_hats,0.975),3),']\n',sep='')
  }

if ( do_sim_tests )
  {
  cat('Simulations:\n')
  S_T_hats = unlist(lapply(all_sim_stats, function(x) x$S_T_hat))
  E_T_hats = unlist(lapply(all_sim_stats, function(x) x$E_T_hat))
  X_T_hats = unlist(lapply(all_sim_stats, function(x) x$X_T_hat))
  p_T_hats = unlist(lapply(all_sim_stats, function(x) x$p_T_hat))
  cat('    S_T_hat = ',round(mean(S_T_hats),1),' [',round(quantile(S_T_hats,0.025),1),', ',round(quantile(S_T_hats,0.975),1),']\n',sep='')
  cat('    E_T_hat = ',round(mean(E_T_hats),1),' [',round(quantile(E_T_hats,0.025),1),', ',round(quantile(E_T_hats,0.975),1),']\n',sep='')
  cat('    X_T_hat = ',round(mean(X_T_hats),1),' [',round(quantile(X_T_hats,0.025),1),', ',round(quantile(X_T_hats,0.975),1),']\n',sep='')
  cat('    p_T_hat = ',round(mean(p_T_hats),3),' [',round(quantile(p_T_hats,0.025),3),', ',round(quantile(p_T_hats,0.975),3),']\n',sep='')
  }
}
