library(tidyverse)
library(scales)

# N=100
# df = data.frame(
#   x=c(seq(0,1,0.1), 1:N)
# ) %>% 
#   mutate(
#     y=2.05^(x/20+rnorm(N+11)) * 50 * rnorm(N+11)^3
#   )


## basic
lal_trans_transform <- function(x) case_when(
  x < -1 ~ -log10(abs(x)) - 1,
  x > 1 ~ log10(x) + 1,
  TRUE ~ x
)

lal_trans_inverse <- function(x) case_when(
  x < -1 ~ -10^(abs(x+1)),
  x > 1 ~ 10^(x-1),
  TRUE ~ x
)

lal_trans = trans_new(
  'lal',
  transform = lal_trans_transform,
  inverse = lal_trans_inverse,
  breaks = function(x) {
    x = x[is.finite(x)]
    
    rng = range(x)
    if (rng[1] < -1){
      min_val = -ceiling(log10(abs(rng[1])+1)) - 1
    } else if (rng[1] < 0){
      min_val = -1
    } else if (rng[1] < 1){
      min_val = 0
    } else {
      min_val = ceiling(log10(rng[1])-1) - 1
    }
    
    if (rng[2] > 1){
      max_val = floor(log10(abs(rng[2]) + 1)) + 1
    } else if (rng[2] > 0){
      max_val = 1
    } else if (rng[2] > -1){
      max_val = 0
    } else {
      max_val = -floor(log10(abs(rng[1]))-1) + 1
    }

    breaks = lal_trans_inverse(as.numeric(seq.int(min_val, max_val)))
    return(breaks)
  }
)

#ggplot(df) + 
#  geom_point(aes(x=x,y=y)) + 
#  scale_y_continuous(trans=lal_trans)


## custom

## make_lal_trans
# Makes a log/absolute log trans object
# where values > threshold use log scale, 
#       values < -threshold use -log scale
#       and values between use linear scale
# @param name - name to use for scale
# @param threshold - threshold magnitude for linear values
# @param exponent - exponent to use for log scale
# @param threshold_scale - if provided, will give the linear region on either side of 0 this much weight 
#                          vs. a unit change in the exponent
# @param return_func_list - return transform and inverse transform functions, in addition to trans object and input parameters
# @param force_thresholds_in_breaks - force the threshold values to be included in the breaks
make_lal_trans <- function(
  name, 
  threshold=1, 
  exponent=10, 
  threshold_scale = NA, 
  return_func_list=FALSE,
  max_breaks=15,
  force_thresholds_in_breaks=FALSE
){
  logf <- function(x) log(x)/log(exponent)
  expf <- function(x) exponent ^ x
  
  if (is.na(threshold_scale)){
    threshold_offset = 0
    threshold_multiplier = 1
  } else {
    threshold_offset = threshold_scale-threshold
    threshold_multiplier = threshold_scale/threshold
  }
  
  cust_lal_transform <- function(x){
    case_when(
      x < -threshold ~ -logf(abs(x)) + logf(threshold) - threshold - threshold_offset,
      x > threshold ~ logf(x) - logf(threshold) + threshold + threshold_offset,
      TRUE ~ x * threshold_multiplier
    )
    
  }
  
  cust_lal_inverse <- function(x){
    case_when(
      x < -threshold * threshold_multiplier ~ -expf(abs(x) - threshold + logf(threshold) - threshold_offset),
      x > threshold * threshold_multiplier ~ expf(x - threshold + logf(threshold) - threshold_offset),
      TRUE ~ x/threshold_multiplier
    )
  }
  
  nt = trans_new(
    name,
    transform = cust_lal_transform,
    inverse = cust_lal_inverse,
    breaks = function(x) {
      x = x[is.finite(x)]
      
      rng = range(x)
      if (rng[1] < -threshold){
        min_val = -ceiling(logf(abs(rng[1])+1)) - 1
      } else if (rng[1] < 0){
        min_val = -threshold
      } else if (rng[1] < threshold){
        min_val = 0
      } else {
        min_val = ceiling(logf(rng[1])-1) - 1
      }
      
      if (rng[2] > threshold){
        max_val = floor(logf(abs(rng[2]) + 1)) + 1
      } else if (rng[2] > 0){
        max_val = 1
      } else if (rng[2] > -threshold){
        max_val = 0
      } else {
        max_val = -floor(logf(abs(rng[1]))-1) + 1
      }
      
      if (min_val < 0){
        lower_breaks = seq.int(min_val - threshold - threshold_offset, min(0, max_val))
      } else {
        lower_breaks = numeric(0)
      }
      if (max_val > 0){
        upper_breaks = seq.int(max_val + threshold + threshold_offset, max(0, min_val)) %>% as.numeric()
      } else {
        upper_breaks = numeric(0)
      }
      
      breaks = c(lower_breaks, upper_breaks)
      
      if (between(0, min_val, max_val) | any(abs(c(min_val, max_val)) < threshold)){
        breaks = c(breaks, 0)
      }
      
      breaks = sort(unique(breaks))
      
      breaks = cust_lal_inverse(breaks)
      if (length(breaks) > max_breaks){
        n_breaks = length(breaks)
        prop_breaks = max_breaks/n_breaks
        factor = ceiling(1/prop_breaks)
        idx = 1:n_breaks
        if (0 %in% breaks)
          z_idx = which(breaks==0)
        else
          z_idx = 0
        breaks = breaks[idx %% factor == z_idx %% factor]
      }
      if (force_thresholds_in_breaks){
        breaks = sort(c(breaks, -threshold, threshold))
      }
      
      return(breaks)
    }
  )
  if (return_func_list){
    return(list(
      trans=nt,
      transform=cust_lal_transform,
      inverse=cust_lal_inverse,
      name=name,
      threshold=threshold,
      exponent=exponent,
      threshold_scale=threshold_scale
    ))
  }
  return(nt)
}

exp_labeller <- function(exponent, digits=0){
  function(x)
  case_when(
    x == 0 ~ '0',
    x < 0 ~ sprintf('-%s^%s', exponent, round(log(abs(x)), digits=digits)),
    x > 0 ~ sprintf('%s^%s', exponent, round(log(x), digits=digits))
  )
}


#my_new_trans = make_lal_trans(
#  'trexp',
#  1/3,
#  3,
#  threshold_scale=1,
#  TRUE
#)

#my_new_trans[['inverse']](my_new_trans[['transform']](as.numeric(-10:10)))
#my_new_trans[['inverse']](my_new_trans[['transform']](seq(-3.5,1, 0.125)))

#my_new_trans[['transform']](-1/2)
#my_new_trans[['transform']](1/2)



#ggplot(df) + 
#  geom_point(aes(x=x,y=as.numeric(x))) + 
#  geom_point(aes(x=x,y=-as.numeric(x))) + 
#  scale_y_continuous(trans=my_new_trans[['trans']], label=exp_labeller(3))
