## multivariate rank histograms

# check fitting dimensions of forecast and observation objects for computing multivariate rank histograms
checkInput_rhfunctions <- function(fc, obs){
  dimfc <- dim(fc)
  dimobs <- dim(obs)
  if(dimfc[1] != dimobs[1] | dimfc[3] != dimobs[2]){
    stop("dimensions of 'fc' and 'obs' do not match")
  }
}

# Average rank histogram, proposed by Thorarinsdottir et al. (2016)
avg_rh <- function(fc, obs){
  checkInput_rhfunctions(fc, obs)
  n <- dim(fc)[1]
  m <- dim(fc)[2]
  d <- dim(fc)[3]
  ranks <- vector(length = n)
  for(nn in 1:n){
    univ_ranks <- apply(cbind(obs[nn, ], t(fc[nn,,])), 1, rank) 
    preranks <- apply(univ_ranks, 1, mean)
    ranks[nn] <- rank(preranks, ties="random")[1]
  }
  return(ranks)
}

# Band depth rank histogram, proposed by Thorarinsdottir et al. (2016)

bd.rank.gen <- function(x){
  d <- dim(x)
  x.prerank <- array(data = NA, dim = d)
  for(i in 1:d[1]){
    z <- x[i, ]
    tmp.ranks <- rank(z)
    x.prerank[i, ] <- tmp.ranks * {d[1] - tmp.ranks} + {tmp.ranks - 1} *
      sapply(z, function(x, z) sum(x == z), z = z)
  }
  x.rank <- apply(x.prerank, 2, mean)
  x.rank <- rank(x.rank, ties = "random")
  return(x.rank)
} 

bd_rh <- function(fc, obs){
  checkInput_rhfunctions(fc, obs)
  n <- dim(fc)[1]
  m <- dim(fc)[2]
  d <- dim(fc)[3]
  ranks <- vector(length = n)
  for(nn in 1:n){    
    ranks[nn] <- bd.rank.gen(cbind(obs[nn, ], t(fc[nn,,])))[1]
  }
  return(ranks)
}

