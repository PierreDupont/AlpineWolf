mcmc.CatSPIM_flat<-function (data, niter = 2400, nburn = 1200, nthin = 5, M = 200, 
          inits = NA, obstype = "poisson", nswap = NA, proppars = list(lam0 = 0.05, 
                                                                       sigma = 0.1, sx = 0.2, sy = 0.2), keepACs = FALSE, keepGamma = FALSE, 
          keepG = FALSE, IDup = "Gibbs", priors = NA, tf = NA) 
{
  library(abind)
  y.obs <- data$y.obs
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- dim(y.obs)[3]
  ncat = data$IDlist$ncat
  nallele = data$IDlist$nallele
  IDcovs = data$IDlist$IDcovs
  n.samples = sum(y.obs)
  constraints = data$constraints
  G.obs = data$G.obs
  if (!is.matrix(G.obs)) {
    G.obs = matrix(G.obs)
  }
  if (!is.list(IDcovs)) {
    stop("IDcovs must be a list")
  }
  nlevels = unlist(lapply(IDcovs, length))
  if (ncol(G.obs) != ncat) {
    stop("G.obs needs ncat number of columns")
  }
  if (!all(is.na(priors))) {
    if (!is.list(priors)) 
      stop("priors must be a list")
    if (!all(names(priors) == c("sigma"))) 
      stop("priors list element name must be sigma")
    warning("Using Uniform prior for sigma")
    usePriors = TRUE
  }
  else {
    warning("No sigma prior entered, using uniform(0,infty).")
    usePriors = FALSE
  }
  if (length(dim(y.obs)) != 3) {
    stop("dim(y.obs) must be 3. Reduced to 2 during initialization")
  }
  if (is.na(nswap)) {
    nswap = n.samples/2
    warning("nswap not specified, using n.samples/2")
  }
  if (!IDup %in% c("MH", "Gibbs")) {
    stop("IDup must be MH or Gibbs")
  }
  if (IDup == "MH") {
  }
  if (obstype == "bernoulli" & IDup == "Gibbs") {
    stop("Must use MH IDup for bernoulli data")
  }
  if ("vertices" %in% names(data)) {
    vertices = data$vertices
    useverts = TRUE
    xlim = c(min(vertices[, 1]), max(vertices[, 1]))
    ylim = c(min(vertices[, 2]), max(vertices[, 2]))
  }
  else if ("buff" %in% names(data)) {
    buff <- data$buff
    xlim <- c(min(X[, 1]), max(X[, 1])) + c(-buff, buff)
    ylim <- c(min(X[, 2]), max(X[, 2])) + c(-buff, buff)
    vertices = cbind(xlim, ylim)
    useverts = FALSE
  }
  else {
    stop("user must supply either 'buff' or 'vertices' in data object")
  }
  if (!any(is.na(tf))) {
    if (any(tf > K)) {
      stop("Some entries in tf are greater than K.")
    }
    if (is.null(dim(tf))) {
      if (length(tf) != J) {
        stop("tf vector must be of length J if summing k dimension over traps.")
      }
      K2D = matrix(rep(tf, M), nrow = M, ncol = J, byrow = TRUE)
    }
    else {
      K2D = rowSums(tf)
      K2D = matrix(rep(K2D, M), nrow = M, ncol = J, byrow = TRUE)
    }
  }
  else {
    tf = rep(K, J)
    K2D = matrix(rep(tf, M), nrow = M, ncol = J, byrow = TRUE)
  }
  psi <- inits$psi
  lam0 <- inits$lam0
  sigma <- inits$sigma
  gamma = inits$gamma
  if (!is.list(gamma)) {
    stop("inits$gamma must be a list")
  }
  if (is.null(constraints)) {
    constraints = matrix(1, nrow = n.samples, ncol = n.samples)
    for (i in 1:n.samples) {
      for (j in 1:n.samples) {
        guys1 = which(G.obs[i, ] != 0)
        guys2 = which(G.obs[j, ] != 0)
        comp = guys1[which(guys1 %in% guys2)]
        if (any(G.obs[i, comp] != G.obs[j, comp])) {
          constraints[i, j] = 0
        }
      }
    }
  }
  if (nrow(constraints) != ncol(constraints)) {
    stop("identity constraint matrix needs to be symmetric")
  }
  binconstraints = FALSE
  if (obstype == "bernoulli") {
    idx = t(apply(y.obs, 1, function(x) {
      which(x > 0, arr.ind = TRUE)
    }))
    for (i in 1:n.samples) {
      for (j in 1:n.samples) {
        if (i != j) {
          if (all(idx[i, 1:2] == idx[j, 1:2])) {
            constraints[i, j] = 0
            constraints[j, i] = 0
            binconstraints = TRUE
          }
        }
      }
    }
  }
  y.true = array(0, dim = c(M, J, K))
  ID = rep(NA, n.samples)
  idx = 1
  for (i in 1:n.samples) {
    if (idx > M) {
      stop("Need to raise M to initialize y.true")
    }
    if (K > 1) {
      traps = which(rowSums(y.obs[i, , ]) > 0)
    }
    else {
      traps = which(y.obs[i, , ] > 0)
    }
    y.true2D = apply(y.true, c(1, 2), sum)
    if (length(traps) == 1) {
      cand = which(y.true2D[, traps] > 0)
    }
    else {
      cand = which(rowSums(y.true2D[, traps]) > 0)
    }
    if (length(cand) > 0) {
      if (length(cand) > 1) {
        cand = cand[1]
      }
      cands = which(ID %in% cand)
      if (all(constraints[i, cands] == 1)) {
        y.true[cand, , ] = y.true[cand, , ] + y.obs[i, 
                                                    , ]
        ID[i] = cand
      }
      else {
        y.true[idx, , ] = y.obs[i, , ]
        ID[i] = idx
        idx = idx + 1
      }
    }
    else {
      y.true[idx, , ] = y.obs[i, , ]
      ID[i] = idx
      idx = idx + 1
    }
  }
  checkID = unique(ID)
  for (i in 1:length(checkID)) {
    idx = which(ID == checkID[i])
    if (!all(constraints[idx, idx] == 1)) {
      stop("ID initialized improperly")
    }
  }
  y.true2D = apply(y.true, c(1, 2), sum)
  known.vector = c(rep(1, max(ID)), rep(0, M - max(ID)))
  z = 1 * (apply(y.true2D, 1, sum) > 0)
  add = M * (0.5 - sum(z)/M)
  if (add > 0) {
    z[sample(which(z == 0), add)] = 1
  }
  s <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], 
                                               ylim[2]))
  idx = which(rowSums(y.true) > 0)
  for (i in idx) {
    trps <- matrix(X[y.true2D[i, ] > 0, 1:2], ncol = 2, byrow = FALSE)
    if (nrow(trps) > 1) {
      s[i, ] <- c(mean(trps[, 1]), mean(trps[, 2]))
    }
    else {
      s[i, ] <- trps
    }
  }
  if (useverts == TRUE) {
    inside = rep(NA, nrow(s))
    for (i in 1:nrow(s)) {
      inside[i] = inout(s[i, ], vertices)
    }
    idx = which(inside == FALSE)
    if (length(idx) > 0) {
      for (i in 1:length(idx)) {
        while (inside[idx[i]] == FALSE) {
          s[idx[i], ] = c(runif(1, xlim[1], xlim[2]), 
                          runif(1, ylim[1], ylim[2]))
          inside[idx[i]] = inout(s[idx[i], ], vertices)
        }
      }
    }
  }
  y.true = y.true2D
  y.obs = apply(y.obs, c(1, 2), sum)
  D = e2dist(s, X)
  lamd <- lam0 * exp(-D * D/(2 * sigma * sigma))
  ll.y = array(0, dim = c(M, J))
  if (obstype == "bernoulli") {
    pd = 1 - exp(-lamd)
    ll.y = dbinom(y.true, K2D, pd * z, log = TRUE)
  }
  else if (obstype == "poisson") {
    ll.y = dpois(y.true, K2D * lamd * z, log = TRUE)
  }
  ll.y.cand = ll.y
  G.true = matrix(0, nrow = M, ncol = ncat)
  for (i in 1:max(ID)) {
    idx = which(ID == i)
    if (length(idx) == 1) {
      G.true[i, ] = G.obs[idx, ]
    }
    else {
      if (ncol(G.obs) > 1) {
        G.true[i, ] = apply(G.obs[idx, ], 2, max)
      }
      else {
        G.true[i, ] = max(G.obs[idx, ])
      }
    }
  }
  G.latent = G.true == 0
  for (j in 1:ncat) {
    fix = G.true[, j] == 0
    G.true[fix, j] = sample(IDcovs[[j]], sum(fix), replace = TRUE, 
                            prob = gamma[[j]])
  }
  nstore = (niter - nburn)/nthin
  if (nburn%%nthin != 0) {
    nstore = nstore + 1
  }
  out <- matrix(NA, nrow = nstore, ncol = 5)
  dimnames(out) <- list(NULL, c("lam0", "sigma", 
                                "N", "n", "psi"))
  if (keepACs) {
    sxout <- syout <- zout <- matrix(NA, nrow = nstore, ncol = M)
    IDout = matrix(NA, nrow = nstore, ncol = length(ID))
  }
  idx = 1
  if (keepGamma) {
    gammaOut = vector("list", ncat)
    for (i in 1:ncat) {
      gammaOut[[i]] = matrix(NA, nrow = nstore, ncol = nlevels[i])
      colnames(gammaOut[[i]]) = paste("Lo", i, "G", 
                                      1:nlevels[i], sep = "")
    }
  }
  if (keepG) {
    Gout = array(NA, dim = c(niter, M, ncat))
  }
  for (iter in 1:niter) {
    if (obstype == "bernoulli") {
      llysum = sum(ll.y)
      lam0.cand <- rnorm(1, lam0, proppars$lam0)
      if (lam0.cand > 0) {
        lamd.cand <- lam0.cand * exp(-D * D/(2 * sigma * 
                                               sigma))
        pd.cand = 1 - exp(-lamd.cand)
        ll.y.cand = dbinom(y.true, K2D, pd.cand * z, 
                           log = TRUE)
        llycandsum = sum(ll.y.cand)
        if (runif(1) < exp(llycandsum - llysum)) {
          lam0 <- lam0.cand
          lamd = lamd.cand
          pd = pd.cand
          ll.y = ll.y.cand
          llysum = llycandsum
        }
      }
      sigma.cand <- rnorm(1, sigma, proppars$sigma)
      if (sigma.cand > 0) {
        lamd.cand <- lam0 * exp(-D * D/(2 * sigma.cand * 
                                          sigma.cand))
        pd.cand = 1 - exp(-lamd.cand)
        ll.y.cand = dbinom(y.true, K2D, pd.cand * z, 
                           log = TRUE)
        llycandsum = sum(ll.y.cand)
        if (usePriors) {
          prior.curr = dunif(sigma, priors$sigma[1], 
                              priors$sigma[2], log = TRUE)
          prior.cand = dunif(sigma.cand, priors$sigma[1], 
                              priors$sigma[2], log = TRUE)
        }
        else {
          prior.curr = prior.cand = 0
        }
        if (runif(1) < exp((llycandsum + prior.cand) - 
                           (llysum + prior.curr))) {
          sigma <- sigma.cand
          lamd = lamd.cand
          pd = pd.cand
          ll.y = ll.y.cand
        }
      }
    }
    else {
      llysum = sum(ll.y)
      lam0.cand <- rnorm(1, lam0, proppars$lam0)
      if (lam0.cand > 0) {
        lamd.cand <- lam0.cand * exp(-D * D/(2 * sigma * 
                                               sigma))
        ll.y.cand = dpois(y.true, K2D * lamd.cand * z, 
                          log = TRUE)
        llycandsum = sum(ll.y.cand)
        if (runif(1) < exp(llycandsum - llysum)) {
          lam0 <- lam0.cand
          lamd = lamd.cand
          ll.y = ll.y.cand
          llysum = llycandsum
        }
      }
      sigma.cand <- rnorm(1, sigma, proppars$sigma)
      if (sigma.cand > 0) {
        lamd.cand <- lam0 * exp(-D * D/(2 * sigma.cand * 
                                          sigma.cand))
        ll.y.cand = dpois(y.true, K2D * lamd.cand * z, 
                          log = TRUE)
        llycandsum = sum(ll.y.cand)
        if (usePriors) {
          prior.curr = dgamma(sigma, priors$sigma[1], 
                              priors$sigma[2], log = TRUE)
          prior.cand = dgamma(sigma.cand, priors$sigma[1], 
                              priors$sigma[2], log = TRUE)
        }
        else {
          prior.curr = prior.cand = 0
        }
        if (runif(1) < exp((llycandsum + prior.cand) - 
                           (llysum + prior.curr))) {
          sigma <- sigma.cand
          lamd = lamd.cand
          ll.y = ll.y.cand
        }
      }
    }
    if (IDup == "Gibbs") {
      up = sample(1:n.samples, nswap, replace = FALSE)
      for (l in up) {
        nj = which(y.obs[l, ] > 0)
        idx2 = which(G.obs[l, ] != 0)
        if (length(idx2) > 1) {
          possible = which(z == 1 & apply(G.true[, idx2], 
                                          1, function(x) {
                                            all(x == G.obs[l, idx2])
                                          }))
        }
        else if (length(idx2) == 1) {
          possible = which(z == 1 & G.true[, idx2] == 
                             G.obs[l, idx2])
        }
        else {
          possible = which(z == 1)
        }
        njprobs = lamd[, nj]
        njprobs[setdiff(1:M, possible)] = 0
        njprobs = njprobs/sum(njprobs)
        newID = sample(1:M, 1, prob = njprobs)
        if (ID[l] != newID) {
          swapped = c(ID[l], newID)
          y.true[ID[l], ] = y.true[ID[l], ] - y.obs[l, 
                                                    ]
          y.true[newID, ] = y.true[newID, ] + y.obs[l, 
                                                    ]
          ID[l] = newID
          ll.y[swapped, ] = dpois(y.true[swapped, ], 
                                  K2D[swapped, ] * lamd[swapped, ], log = TRUE)
        }
      }
    }
    else {
      up = sample(1:n.samples, nswap, replace = FALSE)
      y.cand = y.true
      for (l in up) {
        nj = which(y.obs[l, ] > 0)
        idx2 = which(G.obs[l, ] != 0)
        if (length(idx2) > 1) {
          possible = which(z == 1 & apply(G.true[, idx2], 
                                          1, function(x) {
                                            all(x == G.obs[l, idx2])
                                          }))
        }
        else if (length(idx2) == 1) {
          possible = which(z == 1 & G.true[, idx2] == 
                             G.obs[l, idx2])
        }
        else {
          possible = which(z == 1)
        }
        if (binconstraints) {
          legal = rep(TRUE, length(possible))
          for (i in 1:length(possible)) {
            check = which(ID == possible[i])
            if (length(check) > 0) {
              if (any(constraints[l, check] == 0)) {
                legal[i] = FALSE
              }
            }
          }
          possible = possible[legal]
        }
        njprobs = lamd[, nj]
        njprobs[setdiff(1:M, possible)] = 0
        njprobs = njprobs/sum(njprobs)
        newID = ID
        newID[l] = sample(1:M, 1, prob = njprobs)
        if (ID[l] == newID[l]) 
          next
        swapped = c(ID[l], newID[l])
        propprob = njprobs[swapped[2]]
        backprob = njprobs[swapped[1]]
        y.cand[ID[l], ] = y.true[ID[l], ] - y.obs[l, 
                                                  ]
        y.cand[newID[l], ] = y.true[newID[l], ] + y.obs[l, 
                                                        ]
        focalprob = (sum(ID == ID[l])/n.samples) * (y.true[ID[l], 
                                                           nj]/sum(y.true[ID[l], ]))
        focalbackprob = (sum(newID == newID[l])/n.samples) * 
          (y.cand[newID[l], nj]/sum(y.cand[newID[l], 
                                           ]))
        if (obstype == "poisson") {
          ll.y.cand[swapped, ] = dpois(y.cand[swapped, 
                                              ], K2D[swapped, ] * lamd[swapped, ], log = TRUE)
        }
        else {
          ll.y.cand[swapped, ] = dbinom(y.cand[swapped, 
                                               ], K2D[swapped, ], pd[swapped, ], log = TRUE)
        }
        if (runif(1) < exp(sum(ll.y.cand[swapped, ]) - 
                           sum(ll.y[swapped, ])) * (backprob/propprob) * 
            (focalbackprob/focalprob)) {
          y.true[swapped, ] = y.cand[swapped, ]
          ll.y[swapped, ] = ll.y.cand[swapped, ]
          ID[l] = newID[l]
        }
      }
    }
    known.vector = 1 * (rowSums(y.true) > 0)
    G.true.tmp = matrix(0, nrow = M, ncol = ncat)
    for (i in unique(ID)) {
      idx2 = which(ID == i)
      if (length(idx2) == 1) {
        G.true.tmp[i, ] = G.obs[idx2, ]
      }
      else {
        if (ncol(G.obs) > 1) {
          G.true.tmp[i, ] = apply(G.obs[idx2, ], 2, max)
        }
        else {
          G.true.tmp[i, ] = max(G.obs[idx2, ])
        }
      }
    }
    G.latent = G.true.tmp == 0
    for (j in 1:ncat) {
      swap = G.latent[, j]
      G.true[swap, j] = sample(IDcovs[[j]], sum(swap), 
                               replace = TRUE, prob = gamma[[j]])
    }
    for (j in 1:ncat) {
      x = rep(NA, nlevels[[j]])
      for (k in 1:nlevels[[j]]) {
        x[k] = sum(G.true[z == 1, j] == k)
      }
      gam = rgamma(rep(1, nlevels[[j]]), 1 + x)
      gamma[[j]] = gam/sum(gam)
    }
    if (obstype == "poisson") {
      pd = 1 - exp(-lamd)
    }
    pbar = (1 - pd)^K2D
    prob0 <- exp(rowSums(log(pbar)))
    fc <- prob0 * psi/(prob0 * psi + 1 - psi)
    z[known.vector == 0] <- rbinom(sum(known.vector == 0), 
                                   1, fc[known.vector == 0])
    if (obstype == "bernoulli") {
      ll.y = dbinom(y.true, K2D, pd * z, log = TRUE)
    }
    else {
      ll.y = dpois(y.true, K2D * lamd * z, log = TRUE)
    }
    psi = rbeta(1, 1 + sum(z), 1 + M - sum(z))
    for (i in 1:M) {
      Scand <- c(rnorm(1, s[i, 1], proppars$sx), rnorm(1, 
                                                       s[i, 2], proppars$sy))
      if (useverts == FALSE) {
        inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & 
          Scand[2] < ylim[2] & Scand[2] > ylim[1]
      }
      else {
        inbox = inout(Scand, vertices)
      }
      if (inbox) {
        dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - 
                                                X[, 2])^2)
        lamd.cand[i, ] <- lam0 * exp(-dtmp * dtmp/(2 * 
                                                     sigma * sigma))
        if (obstype == "bernoulli") {
          pd.cand[i, ] = 1 - exp(-lamd.cand[i, ])
          ll.y.cand[i, ] = dbinom(y.true[i, ], K2D[i, 
                                                   ], pd.cand[i, ] * z[i], log = TRUE)
          if (runif(1) < exp(sum(ll.y.cand[i, ]) - sum(ll.y[i, 
                                                            ]))) {
            s[i, ] = Scand
            D[i, ] = dtmp
            lamd[i, ] = lamd.cand[i, ]
            pd[i, ] = pd.cand[i, ]
            ll.y[i, ] = ll.y.cand[i, ]
          }
        }
        else {
          ll.y.cand[i, ] = dpois(y.true[i, ], K2D[i, 
                                                  ] * lamd.cand[i, ] * z[i], log = TRUE)
          if (runif(1) < exp(sum(ll.y.cand[i, ]) - sum(ll.y[i, 
                                                            ]))) {
            s[i, ] = Scand
            D[i, ] = dtmp
            lamd[i, ] = lamd.cand[i, ]
            ll.y[i, ] = ll.y.cand[i, ]
          }
        }
      }
    }
    if (iter > nburn & iter%%nthin == 0) {
      if (keepACs) {
        sxout[idx, ] <- s[, 1]
        syout[idx, ] <- s[, 2]
        zout[idx, ] <- z
        IDout[idx, ] = ID
      }
      if (keepGamma) {
        for (k in 1:ncat) {
          gammaOut[[k]][idx, ] = gamma[[k]]
        }
      }
      if (keepG) {
        Gout[idx, , ] = G.true
      }
      out[idx, ] <- c(lam0, sigma, sum(z), length(unique(ID)), 
                      psi)
      idx = idx + 1
    }
  }
  if (keepACs & keepGamma) {
    out2 = list(out = out, sxout = sxout, syout = syout, 
                zout = zout, IDout = IDout, gammaOut = gammaOut)
  }
  else if (keepACs & !keepGamma) {
    out2 = list(out = out, sxout = sxout, syout = syout, 
                zout = zout, IDout = IDout)
  }
  else if (!keepACs & keepGamma) {
    out2 = list(out = out, gammaOut = gammaOut)
  }
  else {
    out2 = list(out = out)
  }
  if (keepG) {
    out2[[length(out2) + 1]] = Gout
    names(out2)[length(out2)] = "Gout"
  }
  out2
}