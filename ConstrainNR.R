## x if x>0 or 0 if x<=0, for spline definition
pos <- function(x){
  pos1 <- function(x){
    if(x>0){ x } else { 0 }
  }
  sapply(x, pos1)
}

## Spline function per Andersson TM, Dickman PW, Eloranta S, Lambert PC. Estimating and modelling cure in population-based cancer studies within the framework of flexible parametric survival models. 
RevResCubicSpline <- function(x, k=NULL, nk, int=F, pre=NULL, drop=F){ #nk inlcudes external knots
  if(is.null(k)){
    k <- quantile(x, seq(0,1,length.out=nk))
  }
  if(!(is.null(pre))){
    pre <- paste0(pre,".")
  }
  
  K <- length(k)
  
  s <- data.frame(s1=x)
  
  if(int){s <- cbind(1,s)}
  
  if(K>1){
    for(j in 1:(K-2)){
      l <- (k[K-j] - k[1])/(k[K]-k[1])
      
      s[[paste0("s",j)]] <- pos( (k[K-j]-x)^3 ) - l*pos( (k[K]-x)^3 ) - (1-l)*pos( (k[1]-x)^3 )
    }
    
    s[[paste0("s",K-1)]] <- x
  }
  
  s <- as.matrix(s)
  
  colnames(s) <- paste0(pre, colnames(s))
  
  if(drop){
    return(s[,-ncol(s), drop=F])
  }else{
    return(s)
  }
  
}

## NR.beta function adapted from survPen (Mathieu Fauvernier [aut, cre], Laurent Roche [aut], Laurent Remontet [aut], Zoe Uhry [ctb], Nadine Bossard [ctb])
## changes identified with '## ADDED: ####'
NR.beta_cons <- function (build, beta.ini, detail.beta, max.it.beta = 200, tol.beta = 1e-04, cons) {
  type <- build$type
  X <- build$X
  X.GL <- build$X.GL
  event <- build$event
  expected <- build$expected
  leg <- build$leg
  n.legendre <- build$n.legendre
  t1 <- build$t1
  t0 <- build$t0
  tm <- build$tm
  S <- build$S
  p <- build$p
  is.pwcst <- build$is.pwcst
  pwcst.weights <- build$pwcst.weights
  k = 1
  ll.pen = 100
  ll.pen.old = 1
  if (length(beta.ini) == 1) 
    beta.ini <- rep(beta.ini, p)
  if (length(beta.ini) != p) 
    stop("message NR.beta_cons: the length of beta.ini does not equal the number of regression parameters")
  betaold <- beta.ini
  beta1 <- betaold
  if (detail.beta) {
    cat("---------------------------------------------------------------------------------------", 
        "\n", "Beginning regression parameter estimation", 
        "\n", "\n")
  }
  while (abs(ll.pen - ll.pen.old) > tol.beta | any(abs((beta1 - 
                                                        betaold)/betaold) > tol.beta)) {
    if (k > max.it.beta) {
      stop("message NR.beta_cons: Ran out of iterations (", 
           k, "), and did not converge ")
    }
    if (k >= 2) {
      ll.pen.old <- ll.pen
      betaold <- beta1
    }
    predold = X %vec% betaold
    ftold = exp(predold)
    haz.GL.old <- lapply(1:n.legendre, function(i) exp(X.GL[[i]] %vec% 
                                                         betaold))
    if (is.pwcst) {
      deriv.list <- lapply(1:n.legendre, function(i) X.GL[[i]] * 
                             haz.GL.old[[i]] * pwcst.weights[, i])
    }
    else {
      deriv.list <- lapply(1:n.legendre, function(i) X.GL[[i]] * 
                             haz.GL.old[[i]] * leg$weights[i] * tm)
    }
    f.first <- Reduce("+", deriv.list)
    if (type == "net") {
      grad.unpen.beta <- colSums2(-f.first + (event * 
                                                X * ftold)/(ftold + expected))
    }
    else {
      grad.unpen.beta <- colSums2(-f.first + event * X)
    }
    grad <- grad.unpen.beta - S %vec% betaold
    ## ADDED: ####
    if(is.logical(cons)){
      consTF <- cons
    }else{
      consTF <- rep(F, ncol(X));
      consTF[cons] <- T
    }
    grad[consTF] <- 0
    ##############
    
    deriv.2.list <- lapply(1:n.legendre, function(i) X.GL[[i]] %cross% 
                             (deriv.list[[i]]))
    f.second <- Reduce("+", deriv.2.list)
    if (type == "net") {
      Hess.unpen <- -f.second + X %cross% (event * X * 
                                             expected * ftold/(ftold + expected)^2)
    }
    else {
      Hess.unpen <- -f.second
    }
    Hess <- Hess.unpen - S
    ## ADDED: ####
    Hess[consTF,] <- 0
    Hess[,consTF] <- 0
    ##############
    neg.Hess <- -Hess
    R <- try(chol(neg.Hess), silent = TRUE)
    if (inherits(R, "try-error")) {
      u = 0.001
      cpt.while <- 0
      while (inherits(R, "try-error")) {
        if (cpt.while > 100) {
          stop("message NR.beta_cons: did not succeed in inverting Hessian at iteration ", 
               k)
        }
        R <- try(chol(neg.Hess + u * diag(p)), silent = TRUE)
        u <- 5 * u
        cpt.while <- cpt.while + 1
      }
      if (detail.beta) {
        cat("beta Hessian perturbation, ", cpt.while, 
            "iterations", "\n", "\n")
      }
    }
    Vp <- chol2inv(R)
    if (is.pwcst) {
      integral <- lapply(1:n.legendre, function(i) haz.GL.old[[i]] * 
                           pwcst.weights[, i])
      integral <- Reduce("+", integral)
    }
    else {
      integral <- lapply(1:n.legendre, function(i) haz.GL.old[[i]] * 
                           leg$weights[i])
      integral <- tm * Reduce("+", integral)
    }
    if (type == "net") {
      ll.unpenold <- sum(-integral + event * log(ftold + 
                                                   expected))
    }
    else {
      ll.unpenold <- sum(-integral + event * predold)
    }
    ll.pen.old <- ll.unpenold - 0.5 * sum(betaold * (S %vec% 
                                                       betaold))
    if (is.nan(ll.pen.old)) 
      stop("message NR.beta_cons: convergence issues, cannot evaluate log-likelihood")
    pas <- Vp %vec% grad
    beta1 <- betaold + pas
    pred1 <- X %vec% beta1
    ft1 = exp(pred1)
    haz.GL <- lapply(1:n.legendre, function(i) exp(X.GL[[i]] %vec% 
                                                     beta1))
    if (is.pwcst) {
      integral <- lapply(1:n.legendre, function(i) haz.GL[[i]] * 
                           pwcst.weights[, i])
      integral <- Reduce("+", integral)
    }
    else {
      integral <- lapply(1:n.legendre, function(i) haz.GL[[i]] * 
                           leg$weights[i])
      integral <- tm * Reduce("+", integral)
    }
    if (type == "net") {
      ll.unpen <- sum(-integral + event * log(ft1 + expected))
    }
    else {
      ll.unpen <- sum(-integral + event * pred1)
    }
    ll.pen <- ll.unpen - 0.5 * sum(beta1 * (S %vec% beta1))
    if (is.nan(ll.pen)) {
      ll.pen <- ll.pen.old - 1
    }
    if (ll.pen < ll.pen.old - 0.001) {
      cpt.beta <- 1
      while (ll.pen < ll.pen.old - 0.001) {
        if (cpt.beta > 52) 
          stop("message NR.beta_cons: step has been divided by two 52 times in a row, Log-likelihood could not be optimized")
        cpt.beta <- cpt.beta + 1
        pas <- 0.5 * pas
        beta1 <- betaold + pas
        pred1 <- X %vec% beta1
        ft1 = exp(pred1)
        haz.GL <- lapply(1:n.legendre, function(i) exp(X.GL[[i]] %vec% 
                                                         beta1))
        if (is.pwcst) {
          integral <- lapply(1:n.legendre, function(i) haz.GL[[i]] * 
                               pwcst.weights[, i])
          integral <- Reduce("+", integral)
        }
        else {
          integral <- lapply(1:n.legendre, function(i) haz.GL[[i]] * 
                               leg$weights[i])
          integral <- tm * Reduce("+", integral)
        }
        if (type == "net") {
          ll.unpen <- sum(-integral + event * log(ft1 + 
                                                    expected))
        }
        else {
          ll.unpen <- sum(-integral + event * pred1)
        }
        ll.pen <- ll.unpen - 0.5 * sum(beta1 * (S %vec% 
                                                  beta1))
        if (is.nan(ll.pen)) {
          ll.pen <- ll.pen.old - 1
        }
      }
    }
    if (detail.beta) {
      cat("iter beta: ", k, "\n", "betaold= ", round(betaold, 
                                                     4), "\n", "beta= ", round(beta1, 4), "\n", "abs((beta-betaold)/betaold)= ", 
          round(abs((beta1 - betaold)/betaold), 5), "\n", 
          "ll.pen.old= ", round(ll.pen.old, 4), "\n", 
          "ll.pen= ", round(ll.pen, 4), "\n", "ll.pen-ll.pen.old= ", 
          round(ll.pen - ll.pen.old, 5), "\n", "\n")
    }
    k = k + 1
  }
  if (detail.beta) {
    cat("\n", "Beta optimization ok, ", k - 1, "iterations", 
        "\n", "--------------------------------------------------------------------------------------", 
        "\n")
  }
  list(beta = beta1, ll.unpen = ll.unpen, ll.pen = ll.pen, 
       haz.GL = haz.GL, iter.beta = k - 1)
}


## survPen.fit function adapted from survPen (Mathieu Fauvernier [aut, cre], Laurent Roche [aut], Laurent Remontet [aut], Zoe Uhry [ctb], Nadine Bossard [ctb])
## NR.beta replaced by NR.beta_cons
survPen.fit_cons <- function (build, data, formula, max.it.beta = 200, beta.ini = NULL, 
                         detail.beta = FALSE, method = "LAML", tol.beta = 1e-04, cons) {
  formula <- stats::as.formula(formula)
  cl <- build$cl
  n <- build$n
  X <- build$X
  S <- build$S
  S.smf <- build$S.smf
  S.tensor <- build$S.tensor
  S.tint <- build$S.tint
  S.rd <- build$S.rd
  smooth.name.smf <- build$smooth.name.smf
  smooth.name.tensor <- build$smooth.name.tensor
  smooth.name.tint <- build$smooth.name.tint
  smooth.name.rd <- build$smooth.name.rd
  S.scale <- build$S.scale
  lambda <- build$lambda
  rank.S <- build$rank.S
  S.list <- build$S.list
  S.F <- build$S.F
  S.F.list <- build$S.F.list
  U.F <- build$U.F
  p <- build$p
  df.para <- build$df.para
  df.smooth <- build$df.smooth
  list.smf <- build$list.smf
  list.tensor <- build$list.tensor
  list.tint <- build$list.tint
  list.rd <- build$list.rd
  nb.smooth <- build$nb.smooth
  t0 <- build$t0
  t1 <- build$t1
  tm <- build$tm
  event <- build$event
  expected <- build$expected
  type <- build$type
  Z.smf <- build$Z.smf
  Z.tensor <- build$Z.tensor
  Z.tint <- build$Z.tint
  is.pwcst <- build$is.pwcst
  if (is.null(build$optim.rho)) 
    build$optim.rho <- 0
  pwcst.weights <- build$pwcst.weights
  pwcst.breaks <- build$pwcst.breaks
  leg <- build$leg
  n.legendre <- build$n.legendre
  X.GL <- build$X.GL
  if (is.null(beta.ini)) {
    beta.ini = c(log(sum(event)/sum(t1)), rep(0, df.para + 
                                                df.smooth - 1))
  }
  if (any(sapply(list.smf, `[`, "by") != "NULL") | any(sapply(list.tensor, 
                                                              `[`, "by") != "NULL") | any(sapply(list.tint, `[`, "by") != 
                                                                                          "NULL") | build$is.pwcst) {
    beta.ini = rep(0, df.para + df.smooth)
  }
  Algo.optim <- NR.beta_cons(build, beta.ini, detail.beta = detail.beta, 
                        max.it.beta = max.it.beta, tol.beta = tol.beta, cons=cons)
  beta <- Algo.optim$beta
  names(beta) <- colnames(X)
  ll.unpen <- Algo.optim$ll.unpen
  ll.pen <- Algo.optim$ll.pen
  haz.GL <- Algo.optim$haz.GL
  iter.beta <- Algo.optim$iter.beta
  pred1 = X %vec% beta
  ft1 = exp(pred1)
  if (is.pwcst) {
    deriv.list <- lapply(1:n.legendre, function(i) X.GL[[i]] * 
                           haz.GL[[i]] * pwcst.weights[, i])
  }
  else {
    deriv.list <- lapply(1:n.legendre, function(i) X.GL[[i]] * 
                           haz.GL[[i]] * tm * leg$weights[i])
  }
  deriv.2.list <- lapply(1:n.legendre, function(i) X.GL[[i]] %cross% 
                           (deriv.list[[i]]))
  f.first <- Reduce("+", deriv.list)
  if (type == "net") {
    grad.unpen.beta <- colSums2(-f.first + X * event * ft1/(ft1 + 
                                                              expected))
  }
  else {
    grad.unpen.beta <- colSums2(-f.first + X * event)
  }
  grad.beta <- grad.unpen.beta - S %vec% beta
  f.second <- Reduce("+", deriv.2.list)
  if (type == "net") {
    Hess.unpen.beta <- -f.second + X %cross% (X * event * 
                                                expected * ft1/(ft1 + expected)^2)
  }
  else {
    Hess.unpen.beta <- -f.second
  }
  neg.Hess.beta <- -Hess.unpen.beta + S
  R <- try(chol(neg.Hess.beta), silent = TRUE)
  Hess.beta.modif <- FALSE
  if (inherits(R, "try-error")) {
    Hess.beta.modif <- TRUE
    eigen.temp <- eigen(neg.Hess.beta, symmetric = TRUE)
    U.temp <- eigen.temp$vectors
    vp.temp <- eigen.temp$values
    vp.temp[which(vp.temp < 1e-07)] <- 1e-07
    R <- try(chol(U.temp %mult% diag(vp.temp) %mult% t(U.temp)), 
             silent = TRUE)
    warning("beta Hessian was perturbed at convergence")
  }
  neg.Hess.beta <- crossprod(R)
  Vp <- chol2inv(R)
  Ve <- -Vp %mult% Hess.unpen.beta %mult% Vp
  rownames(Ve) <- colnames(Ve) <- rownames(Vp) <- colnames(Vp) <- colnames(X)
  if (nb.smooth != 0) {
    edf <- rowSums(-Hess.unpen.beta * Vp)
    edf1 <- 2 * edf - rowSums(t(Vp %mult% Hess.unpen.beta) * 
                                (Vp %mult% Hess.unpen.beta))
    LCV <- -ll.unpen + sum(edf)
    if (method == "LCV") {
      criterion.val <- LCV
    }
    if (sum(lambda) < .Machine$double.eps) {
      log.abs.S <- 0
      M.p <- 0
      sub.S <- S[1:rank.S, 1:rank.S]
    }
    else {
      M.p <- p - rank.S
      sub.S <- S[1:rank.S, 1:rank.S]
      qr.S <- qr(sub.S)
      log.abs.S <- sum(log(abs(diag(qr.S$qr))))
    }
    log.det.Hess.beta <- as.numeric(2 * determinant(R, logarithm = TRUE)$modulus)
    LAML <- -(ll.pen + 0.5 * log.abs.S - 0.5 * log.det.Hess.beta + 
                0.5 * M.p * log(2 * pi))
    if (method == "LAML") {
      criterion.val <- LAML
    }
    if (build$optim.rho == 1) {
      S.beta <- lapply(1:nb.smooth, function(i) S.list[[i]] %vec% 
                         beta)
      deriv.rho.beta <- matrix(0, nrow = nb.smooth, ncol = p)
      GL.temp <- vector("list", nb.smooth)
      for (i in 1:nb.smooth) {
        deriv.rho.beta[i, ] <- (-Vp) %vec% S.beta[[i]]
        GL.temp[[i]] <- lapply(1:n.legendre, function(j) (X.GL[[j]] %vec% 
                                                            deriv.rho.beta[i, ]) * haz.GL[[j]])
      }
      if (type == "net") {
        temp.deriv3 <- (X * ft1 * (-ft1 + expected)/(ft1 + 
                                                       expected)^3)
        temp.deriv4 <- (X * ft1 * (ft1^2 - 4 * expected * 
                                     ft1 + expected^2)/(ft1 + expected)^4)
      }
      else {
        temp.deriv3 <- temp.deriv4 <- matrix(0)
      }
      if (method == "LCV") {
        mat.temp <- -Vp + Ve
        temp.LAML <- vector("list", 0)
        temp.LAML2 <- vector("list", 0)
        inverse.new.S <- matrix(0)
        minus.eigen.inv.Hess.beta <- 0
        Hess.LCV1 <- matrix(0, nb.smooth, nb.smooth)
      }
      if (method == "LAML") {
        mat.temp <- matrix(0)
        eigen.mat.temp <- 0
        deriv.mat.temp <- vector("list", 0)
        deriv.rho.Ve <- vector("list", 0)
        inverse.new.S <- try(solve.default(sub.S), silent = TRUE)
        if (inherits(inverse.new.S, "try-error")) {
          cat("\n", "LU decomposition failed to invert penalty matrix, trying QR", 
              "\n", "set detail.rho=TRUE for details", 
              "\n")
          inverse.new.S <- try(qr.solve(qr.S))
        }
        if (inherits(inverse.new.S, "try-error")) {
          cat("\n", "LU and QR decompositions failed to invert penalty matrix, trying Cholesky", 
              "\n", "set detail.rho=TRUE for details", 
              "\n")
          inverse.new.S <- chol2inv(chol(sub.S))
        }
        temp.LAML <- lapply(1:nb.smooth, function(i) S.list[[i]][1:rank.S, 
                                                                 1:rank.S])
        temp.LAML2 <- lapply(1:nb.smooth, function(i) -inverse.new.S %mult% 
                               temp.LAML[[i]] %mult% inverse.new.S)
      }
      grad.list <- grad_rho(X.GL, GL.temp, haz.GL, deriv.rho.beta, 
                            leg$weights, tm, nb.smooth, p, n.legendre, S.list, 
                            temp.LAML, Vp, S.beta, beta, inverse.new.S, 
                            X, temp.deriv3, event, expected, type, Ve, mat.temp, 
                            method)
      grad.rho <- grad.list$grad_rho
      if (method == "LCV") 
        grad.rho <- grad.rho + deriv.rho.beta %vec% 
        (-grad.unpen.beta)
      deriv.rho.inv.Hess.beta <- grad.list$deriv_rho_inv_Hess_beta
      deriv.rho.Hess.unpen.beta <- grad.list$deriv_rho_Hess_unpen_beta
      deriv2.rho.beta <- lapply(1:nb.smooth, function(i) matrix(0, 
                                                                nrow = nb.smooth, ncol = p))
      for (j in 1:nb.smooth) {
        for (j2 in 1:nb.smooth) {
          deriv2.rho.beta[[j2]][j, ] <- deriv.rho.inv.Hess.beta[[j2]] %vec% 
            S.beta[[j]] - Vp %vec% (S.list[[j]] %vec% 
                                      deriv.rho.beta[j2, ])
          if (j == j2) {
            deriv2.rho.beta[[j2]][j, ] <- deriv2.rho.beta[[j2]][j, 
            ] - Vp %mult% S.beta[[j2]]
          }
        }
      }
      if (method == "LCV") {
        for (j2 in 1:nb.smooth) {
          Hess.LCV1[, j2] <- deriv2.rho.beta[[j2]] %vec% 
            (-grad.unpen.beta) + deriv.rho.beta %vec% 
            (-Hess.unpen.beta %vec% deriv.rho.beta[j2, 
            ])
        }
        deriv.rho.Ve <- lapply(1:nb.smooth, function(j2) -((deriv.rho.inv.Hess.beta[[j2]] %mult% 
                                                              Hess.unpen.beta - Vp %mult% deriv.rho.Hess.unpen.beta[[j2]]) %mult% 
                                                             (-Vp) - Vp %mult% Hess.unpen.beta %mult% deriv.rho.inv.Hess.beta[[j2]]))
        deriv.mat.temp <- lapply(1:nb.smooth, function(j2) deriv.rho.Ve[[j2]] + 
                                   deriv.rho.inv.Hess.beta[[j2]])
        eigen2 <- eigen(mat.temp, symmetric = TRUE)
        eigen.mat.temp <- eigen2$values
      }
      if (method == "LAML") {
        eigen2 <- eigen(Vp, symmetric = TRUE)
        minus.eigen.inv.Hess.beta <- eigen2$values
      }
      Q <- eigen2$vectors
      X.Q <- X %mult% Q
      X.GL.Q <- lapply(1:n.legendre, function(i) X.GL[[i]] %mult% 
                         Q)
      Hess.rho <- Hess_rho(X.GL, X.GL.Q, GL.temp, haz.GL, 
                           deriv2.rho.beta, deriv.rho.beta, leg$weights, 
                           tm, nb.smooth, p, n.legendre, deriv.rho.inv.Hess.beta, 
                           deriv.rho.Hess.unpen.beta, S.list, minus.eigen.inv.Hess.beta, 
                           temp.LAML, temp.LAML2, Vp, S.beta, beta, inverse.new.S, 
                           X, X.Q, temp.deriv3, temp.deriv4, event, expected, 
                           type, Ve, deriv.rho.Ve, mat.temp, deriv.mat.temp, 
                           eigen.mat.temp, method)
      if (method == "LCV") 
        Hess.rho <- Hess.rho + Hess.LCV1
    }
    else {
      grad.rho <- Hess.rho <- deriv.rho.beta <- deriv.rho.inv.Hess.beta <- NULL
    }
  }
  else {
    edf <- edf1 <- p
    LCV <- LAML <- criterion.val <- grad.rho <- Hess.rho <- deriv.rho.beta <- deriv.rho.inv.Hess.beta <- NULL
  }
  optim.rho <- iter.rho <- edf2 <- aic2 <- inv.Hess.rho <- Hess.rho.modif <- Vc.approx <- Vc <- NULL
  res <- list(call = cl, formula = formula, t0.name = build$t0.name, 
              t1.name = build$t1.name, event.name = build$event.name, 
              expected.name = build$expected.name, haz = ft1, coefficients = beta, 
              type = type, df.para = df.para, df.smooth = df.smooth, 
              p = p, edf = edf, edf1 = edf1, edf2 = edf2, aic = 2 * 
                sum(edf) - 2 * ll.unpen, aic2 = aic2, iter.beta = iter.beta, 
              X = X, S = S, S.scale = S.scale, S.list = S.list, S.smf = S.smf, 
              S.tensor = S.tensor, S.tint = S.tint, S.rd = S.rd, smooth.name.smf = smooth.name.smf, 
              smooth.name.tensor = smooth.name.tensor, smooth.name.tint = smooth.name.tint, 
              smooth.name.rd = smooth.name.rd, S.pen = build$S.pen, 
              grad.unpen.beta = grad.unpen.beta, grad.beta = grad.beta, 
              Hess.unpen.beta = Hess.unpen.beta, Hess.beta = -neg.Hess.beta, 
              Hess.beta.modif = Hess.beta.modif, ll.unpen = ll.unpen, 
              ll.pen = ll.pen, deriv.rho.beta = deriv.rho.beta, deriv.rho.inv.Hess.beta = deriv.rho.inv.Hess.beta, 
              lambda = lambda, nb.smooth = nb.smooth, iter.rho = iter.rho, 
              optim.rho = optim.rho, method = method, criterion.val = criterion.val, 
              LCV = LCV, LAML = LAML, grad.rho = grad.rho, Hess.rho = Hess.rho, 
              inv.Hess.rho = inv.Hess.rho, Hess.rho.modif = Hess.rho.modif, 
              Ve = Ve, Vp = Vp, Vc = Vc, Vc.approx = Vc.approx, Z.smf = Z.smf, 
              Z.tensor = Z.tensor, Z.tint = Z.tint, list.smf = list.smf, 
              list.tensor = list.tensor, list.tint = list.tint, list.rd = list.rd, 
              U.F = U.F, is.pwcst = is.pwcst, pwcst.breaks = pwcst.breaks)
  class(res) <- "survPen"
  res
}

## survPen function adapted from survPen package (Mathieu Fauvernier [aut, cre], Laurent Roche [aut], Laurent Remontet [aut], Zoe Uhry [ctb], Nadine Bossard [ctb])
## survPen.fit replaced by survPen.fit_cons
survPen_cons <- function (formula, data, t1, t0 = NULL, event, expected = NULL, 
                     lambda = NULL, rho.ini = NULL, max.it.beta = 200, max.it.rho = 30, 
                     beta.ini = NULL, detail.rho = FALSE, detail.beta = FALSE, 
                     n.legendre = 20, method = "LAML", tol.beta = 1e-04, tol.rho = 1e-04, 
                     step.max = 5, cons=F){
  cl <- match.call()
  if (missing(formula) | missing(data) | missing(t1) | missing(event)) 
    stop("Must have at least a formula, data, t1 and event arguments")
  formula <- stats::as.formula(formula)
  if (!(method %in% c("LAML", "LCV"))) 
    stop("method should be LAML or LCV")
  data <- as.data.frame(unclass(data), stringsAsFactors = TRUE)
  factor.term <- names(data)[sapply(data, is.factor)]
  for (factor.name in names(data)[names(data) %in% factor.term]) {
    data[, factor.name] <- factor(data[, factor.name])
  }
  factor.structure <- lapply(as.data.frame(data[, names(data) %in% 
                                                  factor.term]), attributes)
  names(factor.structure) <- factor.term
  t1.name <- deparse(substitute(t1))
  t1 <- eval(substitute(t1), data)
  t1 <- as.numeric(t1)
  n <- length(t1)
  t0.name <- deparse(substitute(t0))
  t0 <- eval(substitute(t0), data)
  event.name <- deparse(substitute(event))
  event <- eval(substitute(event), data)
  expected.name <- deparse(substitute(expected))
  expected <- eval(substitute(expected), data)
  if (is.null(expected)) {
    type <- "overall"
    expected <- rep(0, n)
  }
  else {
    type <- "net"
    expected <- as.numeric(expected)
  }
  if (is.null(t0)) {
    t0 <- rep(0, n)
  }
  else {
    t0 <- as.numeric(t0)
  }
  if (length(t0) == 1) 
    t0 <- rep(t0, n)
  if (is.null(event)) {
    event <- rep(1, n)
  }
  else {
    event <- as.numeric(event)
  }
  if (any(t0 > t1)) 
    stop("some t0 values are superior to t1 values")
  if (length(t0) != n) 
    stop("t0 and t1 are different lengths")
  if (length(event) != n) 
    stop("event and t1 are different lengths")
  if (length(expected) != n) 
    stop("expected and t1 are different lengths")
  build <- model.cons(formula, lambda, data, t1, t1.name, 
                      t0, t0.name, event, event.name, expected, expected.name, 
                      type, n.legendre, cl, beta.ini)
  if (is.null(lambda)) {
    nb.smooth <- build$nb.smooth
    if (nb.smooth != 0) {
      if (build$is.pwcst) 
        stop("Must not use piecewise constant hazard specification with penalized splines that need smoothing parameter estimation. Please use the cut function")
      if (is.null(rho.ini)) 
        rho.ini <- rep(-1, nb.smooth)
      if (length(rho.ini) != nb.smooth) {
        if (length(rho.ini) == 1) {
          rho.ini <- rep(rho.ini, nb.smooth)
        }
        else {
          stop("number of initial log smoothing parameters incorrect")
        }
      }
      param <- repam(build)
      build <- param$build
      X.ini <- param$X.ini
      S.pen.ini <- param$S.pen.ini
      beta.ini <- build$beta.ini
      build$optim.rho <- 1
      model <- NR.rho(build, rho.ini = rho.ini, data = data, 
                      formula = formula, max.it.beta = max.it.beta, 
                      max.it.rho = max.it.rho, beta.ini = beta.ini, 
                      detail.rho = detail.rho, detail.beta = detail.beta, 
                      nb.smooth = nb.smooth, tol.beta = tol.beta, 
                      tol.rho = tol.rho, step.max = step.max, method = method)
      model <- inv.repam(model, X.ini, S.pen.ini)
      if (model$method == "LAML") 
        model <- cor.var(model)
      model$factor.structure <- factor.structure
      model$converged <- !(model$Hess.beta.modif | model$Hess.rho.modif)
      return(model)
    }
    else {
      build$lambda <- 0
      build$optim.rho <- 0
      model <- survPen.fit_cons(build, data = data, formula = formula, 
                           max.it.beta = max.it.beta, beta.ini = beta.ini, 
                           detail.beta = detail.beta, method = method, 
                           tol.beta = tol.beta, cons=cons)
      model$factor.structure <- factor.structure
      model$converged <- !(model$Hess.beta.modif)
      return(model)
    }
  }
  else {
    param <- repam(build)
    build <- param$build
    X.ini <- param$X.ini
    S.pen.ini <- param$S.pen.ini
    beta.ini <- build$beta.ini
    build$optim.rho <- 0
    model <- survPen.fit_cons(build, data = data, formula = formula, 
                         max.it.beta = max.it.beta, beta.ini = beta.ini, 
                         detail.beta = detail.beta, method = method, tol.beta = tol.beta, cons=cons)
    model <- inv.repam(model, X.ini, S.pen.ini)
    model$factor.structure <- factor.structure
    model$converged <- !(model$Hess.beta.modif)
    return(model)
  }
}

## Perform regression standardisation for survPen object
standSurvPen <- function(object, data, at, times, type=c("haz","surv"), contrast=NULL, n.legendre=30, 
                         conf.int=0.95, overall=F, gran=0.5, expected=NULL, agemap=NULL, diagmap=NULL, othermap=NULL){
  #diagmap list with colname in first entry then diagformat="%Y", scale=365.25 as extra arguments
  
  #map in the main function should map to rates at the event time (ie floor(year+t)), but for standardisation should be bl (ie year - in day format for this one)
  #why difference in date format between the 2 functions?? ie why does stand. need it in day format?
  
  #checks like no haz if overall, no ci if overall, all ages/diags have matching expecetd vals (ie no pts too young or old)
  #others:
  ##things like, overall F if no expecteds, expecteds are rightish form, at names allowed, type argument ok, contrast argument ok, maps all look ok incl all names of diagmap ok (1st in df and rest as allowed), only overall if surv, only 2 at levels if using contrast, cannot have overall T if haz in type
  
  gl_base <- gaussLegendre(n.legendre,-1,1)
  
  if("rmst" %in% type){
    if(length(times)>1){stop("Only 1 time value allowed for rmst calculation.")}
    time <- times 
    
    times <- (time/2)*gl_base$x + time/2
  }
  
  #S*
  sstar <- lapply(seq_along(at), FUN=function(i){
    
    if(object$expected.name=="NULL" & is.null(expected)){#no expected argument
      sstar <- matrix(1, nrow=nrow(data), ncol=length(times))
    } else if(!(overall)){#only looking at excess
      sstar <- matrix(1, nrow=nrow(data), ncol=length(times))
    } else{#overall survival - need to change this if add in overall hazard option???
      
      #UPDATE:::::
      if(is.null(expected)){expected <- cS$expected}
      
      at1 <- at[[i]]
      d <- data
      d[names(at1)] <- rep(unlist(at1), each=nrow(d))
      d$at_id <- i #check col doesn't exist
      
      splits <- unique(c(seq(0,max(data[[object$t1.name]]),gran),max(data[[object$t1.name]])))
      map=c(agemap, diagmap, othermap)
      
      sstar <- sapply(times,FUN=function(time){
        
        df <- d %>% 
          mutate(s_t=time, s_d=T, s_id=1:nrow(data))
        df_sS <- df %>% 
          survSplit(formula=Surv(s_t,s_d)~.,cut=splits)
        
        if(!(is.null(agemap))){df_sS <- df_sS %>% mutate(age_dyn = floor(get(names(agemap))+tstart) )}
        if(!(is.null(diagmap))){
          if(is.null(diagmap[["diagformat"]])){diagformat <- "%Y"}else{diagformat <- diagmap[["diagformat"]]}
          if(is.null(diagmap[["scale"]])){scale <- 365.25}else{scale <- diagmap[["scale"]]}
          df_sS <- df_sS %>% mutate(diag_dyn = format(get(names(diagmap))+tstart*scale, format = diagformat) )
        }
        
        #add warnings here:
        
        if(!(is.null(agemap))){
          df_sS[df_sS[,"age_dyn"] > apply(expected[,unname(unlist(agemap))], MARGIN = 2, max), "age_dyn"] <- apply(expected[,unname(unlist(agemap))], MARGIN = 2, max)
        }
        if(!(is.null(diagmap))){
          df_sS[df_sS[,"diag_dyn"] > apply(expected[,unname(unlist(diagmap))], MARGIN = 2, max), "diag_dyn"] <- apply(expected[,unname(unlist(diagmap))], MARGIN = 2, max)
        }
        
        
        for(j in seq_along(map)){
          m0 <- c(ifelse(is.null(agemap),NULL,"age_dyn"), ifelse(is.null(diagmap),NULL,"diag_dyn"), names(othermap))[j]
          df_sS[df_sS[,m0] > apply(expected[,unname(unlist(map))], MARGIN = 2, max)[j], m0] <- apply(expected[,unname(unlist(map))], MARGIN = 2, max)[j]
        }
        
        df_sS <- merge(x=df_sS, y=expected, by.x=c(ifelse(is.null(agemap),NULL,"age_dyn"), ifelse(is.null(diagmap),NULL,"diag_dyn"), names(othermap)), by.y=unname(unlist(map)), all.x=TRUE) %>% 
          arrange(s_id,tstart) %>% 
          group_by(s_id) %>% 
          mutate(y=s_t - tstart) %>% 
          mutate(ch=cumsum(y*rate)) %>% 
          mutate(surv=exp(-ch))
        
        df_sS[df_sS$s_t==time,]$surv
      })
    }
    colnames(sstar) <- times
    sstar
    
  }) %>% list.rbind()
  
  if(object$expected.name=="NULL" & is.null(expected)){
  } else if(!(overall)){
  } else{
    conf.int <- F 
  }
  
  out <- lapply(times, FUN=function(t0){
    
    ta <- ifelse( isFALSE(conf.int), 0, qnorm(1 - (1-conf.int)/2) )
    
    msurv <- mhaz <- NULL
    
    df <- lapply(seq_along(at), FUN=function(i){
      at1 <- at[[i]]
      d <- data[rep(1:nrow(data), each=n.legendre),]
      d[[object$t1.name]] <- rep((t0/2)*gl_base$x + t0/2, nrow(data))
      d[names(at1)] <- rep(unlist(at1), each=nrow(d))
      d$at_id <- i #check col doesn't exist
      d
    }) %>% list.rbind()
    
    Xmat <- predict(object, type="lpmatrix", newdata=df)
    
    haz <- as.vector( exp(Xmat%*%t(t(object$coefficients))) )
    
    hazmat <- matrix(haz, byrow=T,ncol=n.legendre)
    
    S <- exp(-(t0/2)*colSums(gl_base$w*t(hazmat)))
    
    Gmat <- rowsum( rep(gl_base$w, nrow(data)*length(at))*haz*Xmat , rep(1:(nrow(data)*length(at)), each=n.legendre) )
    
    if("rmst" %in% type){
      msurv <- by(matrix(S*sstar[,as.character(t0)],ncol=1), INDICES = rep(1:length(at), each=nrow(data)), FUN=function(d){
        colMeans(d[,1, drop=F])
      }) %>% list.rbind() %>% as.data.frame()
      colnames(msurv) <- "Est"
      
      if(is.null(names(at))){
        rownames(msurv) <- paste0("at",1:length(at))
      } else {
        rownames(msurv) <- names(at)
      }
      
    }
    
    if("surv" %in% type){
      
      S_G <- cbind(S,Gmat)
      
      G <- by(S_G, INDICES = rep(1:length(at), each=nrow(data)), FUN=function(d){
        colMeans(-(t0/2)*d[,1]*d[,-1])
      }) %>% list.rbind()
      
      var <- G%*%object$Ve%*%t(G)
      
      if(is.null(names(at))){
        colnames(var) <- rownames(var) <- paste0("at",1:length(at))
      } else {
        colnames(var) <- rownames(var) <- names(at)
      }
      
      #msurv
      msurv <- by(matrix(S*sstar[,as.character(t0)],ncol=1), INDICES = rep(1:length(at), each=nrow(data)), FUN=function(d){
        colMeans(d[,1, drop=F])
      }) %>% list.rbind() %>% as.data.frame()
      colnames(msurv) <- "Est"
      rownames(msurv) <- rownames(var)
      
      # var for comp log surv
      msurv$Var <- diag(var)/((msurv$Est^2)*(log(msurv$Est)^2))
      
      msurv$LCI <- exp(-exp( log(-log(msurv$Est)) + ta*sqrt(msurv$Var) ))
      msurv$UCI <- exp(-exp( log(-log(msurv$Est)) - ta*sqrt(msurv$Var) ))
      
      if("diff" %in% contrast){
        #msdiff
        msdiff <- data.frame(Est=msurv$Est[1]-msurv$Est[2])
        rownames(msdiff) <- "diff"
        
        msdiff$Var <- diag(var)[1] + diag(var)[2] -2*var[1,2]
        
        msdiff$LCI <- msdiff$Est - ta*sqrt(msdiff$Var)
        msdiff$UCI <- msdiff$Est + ta*sqrt(msdiff$Var)
        
        msurv <- rbind(msurv, msdiff)
        
      }
      
      if("ratio" %in% contrast){
        #msratio
        msratio <- data.frame(Est=msurv$Est[1]/msurv$Est[2])
        rownames(msratio) <- "ratio"
        
        # var for log ratio
        msratio$Var <- (1/msurv$Est[1]^2)*(diag(var)[1] + var[1,2]) + (1/msurv$Est[2]^2)*(diag(var)[2] + var[1,2])
        
        msratio$LCI <- exp(log(msratio$Est) - ta*sqrt(msratio$Var))
        msratio$UCI <- exp(log(msratio$Est) + ta*sqrt(msratio$Var))
        
        msurv <- rbind(msurv, msratio)
      }
    }
    
    if("haz" %in% type){
      
      df_h <- lapply(seq_along(at), FUN=function(i){
        at1 <- at[[i]]
        d <- data
        d[[object$t1.name]] <- t0
        d[names(at1)] <- rep(unlist(at1), each=nrow(d))
        d$at_id <- i #check col doesn't exist
        d
      }) %>% list.rbind()
      
      Xmat_h <- predict(object, type="lpmatrix", newdata=df_h)
      
      hazs <- as.vector( exp(Xmat_h%*%t(t(object$coefficients))) )
      
      S_h_G_X <- cbind(S, hazs, Gmat, Xmat_h)
      
      G_h <- by(S_h_G_X, INDICES = rep(1:length(at), each=nrow(data)), FUN=function(d){
        (sum(d[,1]) * colSums( -(t0/2)*d[,1]*d[,2]*d[,3:(3+ncol(Gmat)-1)] + d[,1]*d[,2]*d[,(ncol(Gmat)+3):ncol(d)] ) - sum(d[,1]*d[,2]) * colSums( -(t0/2)*d[,1]*d[,3:(3+ncol(Gmat)-1)] )) / sum(d[,1])^2 #(sum(S) * colSums( -(t0/2)*S*hazs*Gmat + S*hazs*Xmat_h ) - sum(S*hazs) * colSums( -(t0/2)*S*Gmat )) / sum(S)^2
      }) %>% list.rbind()
      
      var_h <- G_h%*%object$Ve%*%t(G_h)
      
      if(is.null(names(at))){
        colnames(var_h) <- rownames(var_h) <- paste0("at",1:length(at))
      } else {
        colnames(var_h) <- rownames(var_h) <- names(at)
      }
      
      #mhaz
      mhaz <- by(matrix(c(S,hazs),ncol=2), INDICES = rep(1:length(at), each=nrow(data)), FUN=function(d){
        sum(d[,1]*d[,2])/sum(d[,1])
      }) %>% list.rbind() %>% as.data.frame()
      colnames(mhaz) <- "Est"
      rownames(mhaz) <- rownames(var_h)
      
      # var for log haz
      mhaz$Var <- diag(var_h)/(mhaz$Est^2)
      
      mhaz$LCI <- exp( log(mhaz$Est) - ta*sqrt(mhaz$Var) )
      mhaz$UCI <- exp( log(mhaz$Est) + ta*sqrt(mhaz$Var) )
      
      if("diff" %in% contrast){
        #msdiff
        mhdiff <- data.frame(Est=mhaz$Est[1]-mhaz$Est[2])
        rownames(mhdiff) <- "diff"
        
        mhdiff$Var <- diag(var_h)[1] + diag(var_h)[2] -2*var_h[1,2]
        
        mhdiff$LCI <- mhdiff$Est - ta*sqrt(mhdiff$Var)
        mhdiff$UCI <- mhdiff$Est + ta*sqrt(mhdiff$Var)
        
        mhaz <- rbind(mhaz, mhdiff)
      }
      
      if("ratio" %in% contrast){
        #msratio
        mhratio <- data.frame(Est=mhaz$Est[1]/mhaz$Est[2])
        rownames(mhratio) <- "ratio"
        
        # var for log ratio
        mhratio$Var <- (1/mhaz$Est[1]^2)*(diag(var_h)[1] + var_h[1,2]) + (1/mhaz$Est[2]^2)*(diag(var_h)[2] + var_h[1,2])
        
        mhratio$LCI <- exp(log(mhratio$Est) - ta*sqrt(mhratio$Var))
        mhratio$UCI <- exp(log(mhratio$Est) + ta*sqrt(mhratio$Var))
        
        mhaz <- rbind(mhaz, mhratio)
      }
    }
    
    l <- list(msurv = msurv, mhaz = mhaz)
    
    return(l)
    
  })
  names(out) <- times
  
  if(!(is.null(out[[1]]$msurv))){
    msurv <- lapply(1:nrow(out[[1]]$msurv), FUN=function(i){
      l <- lapply(out, FUN=function(t){
        t$msurv[i,, drop=F]
      }) %>% list.rbind() %>% rownames_to_column()
      colnames(l)[1] <- "t"
      l[,1]<-as.numeric(l[,1])
      l
    }) %>% setNames(rownames(out[[1]]$msurv))
    
    if(!(is.null(msurv$diff))){
      msurv$diff$Contrast <- paste0(names(msurv)[1],"-",names(msurv)[2])
    }
    if(!(is.null(msurv$ratio))){
      msurv$ratio$Contrast <- paste0(names(msurv)[1],"/",names(msurv)[2])
    }
  }else{
    msurv <- NULL
  }
  
  if(!(is.null(out[[1]]$mhaz))){
    mhaz <- lapply(1:nrow(out[[1]]$mhaz), FUN=function(i){
      l <- lapply(out, FUN=function(t){
        t$mhaz[i,]
      }) %>% list.rbind() %>% rownames_to_column()
      colnames(l)[1] <- "t"
      l[,1]<-as.numeric(l[,1])
      l
    }) %>% setNames(rownames(out[[1]]$mhaz))
    
    if(!(is.null(mhaz$diff))){
      mhaz$diff$Contrast <- paste0(names(msurv)[1],"-",names(msurv)[2])
    }
    if(!(is.null(mhaz$ratio))){
      mhaz$ratio$Contrast <- paste0(names(msurv)[1],"/",names(msurv)[2])
    }
  }else{
    mhaz <- NULL
  }
  
  if("rmst" %in% type){
    mrmst <- lapply(msurv, FUN=function(at1){
      Est <- (time/2)*colSums(gl_base$w*t(matrix(at1$Est, byrow=T,ncol=n.legendre)))
      data.frame(t=time,Est=Est)
    })
    lis <- list(mrmst=mrmst)
  }else{
    lis <- list(msurv=msurv, mhaz=mhaz)
  }
  
  return(lis)
  
}
