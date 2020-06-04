dslcd = function (lcd, x, h) {
  nx = length(x)
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  eqx = list()
  eqx$x0 = matrix(0, nrow=nx, ncol=nk)
  for (j in 1:(nk-1)) {
    mu = x + lcd$coef[2, j] * h
    exp = exp(lcd$coef[2, j] * x + lcd$coef[2, j]^2 * h / 2 + lcd$coef[1, j])
    l = (knots[j] - mu) / sqrt(h)
    u = (knots[j+1] - mu) / sqrt(h)
    Phi.diff = pnorm(u) - pnorm(l)
    eqx$x0[, j] = exp * Phi.diff
  }
  smo.f = apply(eqx$x0[, (nk:1)], 1, sum)
  smo.f = smo.f / lcd$C
  return(smo.f)
}

logLik.slcd = function (lcd, x, h) {
  smo.f = dslcd(lcd, x, h)
  return(sum(log(smo.f)))
}

Eqx = function (lcd, x, h) {
  nx = length(x)
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  eqx = list()
  eqx$exp = eqx$x2 = eqx$x1 = eqx$x0 = matrix(0, nrow=nx, ncol=nk)
  for (j in 1:(nk-1)) {
    mu = x + lcd$coef[2, j] * h
    eqx$exp[, j] = exp(lcd$coef[2, j] * x + lcd$coef[2, j]^2 * h / 2 + lcd$coef[1, j])
    if (any(eqx$exp[, j] > 1e+500)) {
      cat("You should use larger h!", "\n")
      break
    }
    l = (knots[j] - mu) / sqrt(h)
    u = (knots[j+1] - mu) / sqrt(h)
    Phi.diff = pnorm(u) - pnorm(l)
    phi.diff = dnorm(l) - dnorm(u)
    chisq.diff = sign(u) * pchisq(u^2, 3) - sign(l) * pchisq(l^2, 3)
    eqx$x0[, j] = eqx$exp[, j] * Phi.diff
    eqx$x1[, j] = eqx$exp[, j] * (mu * Phi.diff + sqrt(h) * phi.diff)
    eqx$x2[, j] = eqx$exp[, j] * (h * chisq.diff / 2 + 2 * sqrt(h) * mu * phi.diff + mu^2 * Phi.diff)
  }
  eqx$x0 = t(apply(eqx$x0[, (nk:1)], 1, cumsum))[, (nk:1)]
  eqx$x1 = t(apply(eqx$x1[, (nk:1)], 1, cumsum))[, (nk:1)]
  eqx$x2 = t(apply(eqx$x2[, (nk:1)], 1, cumsum))[, (nk:1)]
  return(eqx)
}

gradient.delta = function(lcd, x, w, h, delta) {
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  nd = length(delta)
  eqx = Eqx(lcd, x, h)
  cpd = cpx(lcd, delta, order=1)
  efx = lcd$cpk[2, nk - 1] - cpd[2, ] -  delta * (1-cpd[1, ]) # Ef[(X - delta)+]
  dind = lsei::indx(delta, knots)
  dind[dind >= nk] = nk - 1
  grad = numeric(nd)
  for (i in 1:nd) {
    j = dind[i]
    mu = x + lcd$coef[2, j] * h
    l = (delta[i] - mu) / sqrt(h)
    u = (knots[j+1] - mu) / sqrt(h)
    Phi.diff = pnorm(u) - pnorm(l)
    phi.diff = dnorm(l) - dnorm(u)
    eqx0 = eqx$exp[, j] * Phi.diff + eqx$x0[, j+1]
    eqx1 = eqx$exp[, j] * (mu * Phi.diff + sqrt(h) * phi.diff) + eqx$x1[, j+1]
    grad[i] = w %*% ((-eqx1 + delta[i] * eqx0) / eqx$x0[, 1] + efx[i])
  }
  return(grad)
}

maxima.delta = function (lcd, x, w, h, delta, tol=0) {
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  index = match(knots, delta)
  grad = gradient.delta(lcd, x, w, h, delta)
  ii = numeric(0)
  for (i in 1:(nk-1)) {
    grad.between = grad[index[i]:index[i+1]]
    if (max(grad.between) > 0) {
      ima = which.max(grad.between)
      ii = c(ii, index[i] + ima - 1)
    }
  }
  delta = delta[ii]
  g = grad[ii]
  gmax = ifelse(length(g) > 0, max(g), -Inf)
  j = g > tol & ! delta %in% knots
  # maxima.ind = which(diff(diff(g) > 0) == -1) + 1
  # delta = delta[maxima.ind]
  # gmax = g[maxima.ind]
  # j = gmax > tol & ! delta %in% knots
  list(delta=delta[j], gradient=g[j], gmax, gds=grad)
}

gradient.theta = function(lcd, x, w, h) {
  nx = length(x)
  knots = c(lcd$lower, lcd$theta)
  delta = c(0, lcd$theta)
  nk = length(knots)
  # Gradients
  eqx = Eqx(lcd, x, h)
  cpkr = lcd$cpk[, nk] - cbind(0, lcd$cpk[, -nk, drop = FALSE])
  efx = cpkr[2, ] - cpkr[1, ] * delta
  S = matrix(0, nrow=nx, ncol=nk)
  for (j in 1:nk) S[, j] = (delta[j] * eqx$x0[, j] - eqx$x1[, j]) / eqx$x0[, 1] + efx[j]
  grad = drop(t(S) %*% rep(1, nx))
  # Hessian Matrix
  mm = cpkr[3, ] - (delta + rep(delta, rep.int(nk, nk))) * cpkr[2, ] + tcrossprod(delta) * cpkr[1, ]
  mm[upper.tri(mm)] = 0
  mm = mm + t(mm)
  diag(mm) = diag(mm)/2
  H = matrix(0, nrow=nk, ncol=nk)
  for (j in 1:nk) {
    for (k in j:nk) {
      # sign = (-1)^(j != 1) * (-1)^(k != 1)
      eqxl = (eqx$x2[, k] - (delta[j] + delta[k]) * eqx$x1[, k] + delta[j] * delta[k]) * eqx$x0[, 1]
      eqxr = (eqx$x1[, j] - delta[j] * eqx$x0[, j]) * (eqx$x1[, k] - delta[k] * eqx$x0[, k])
      derivative2 = (eqxr - eqxl) / eqx$x0[, 1]^2
      H[j, k] = sum(derivative2)
    }
  }
  H = H + t(H)
  diag(H) = diag(H)/2
  H = H/nx + mm
  return(list(grad = grad, H = H))
}

gradient.beta = function(lcd, x, h) {
  library(miscTools)
  nx = length(x)
  knots = c(lcd$lower, lcd$theta, lcd$upper)
  nk = length(knots)
  grad = list()
  grad$beta0 = grad$beta1 = grad$beta2 = matrix(0, nrow=nx, ncol=nk)
  for (j in 1:(nk-1)) {
    mu = x + lcd$coef[2, j] * h
    exp = exp(lcd$coef[2, j] * x + lcd$coef[2, j]^2 * h / 2 + lcd$coef[1, j])
    l = (knots[j] - mu) / sqrt(h)
    u = (knots[j+1] - mu) / sqrt(h)
    Phi.diff = pnorm(u) - pnorm(l)
    phi.diff = dnorm(l) - dnorm(u)
    dphi.diff = ddnorm(u) - ddnorm(l)
    grad$beta0[, j] = exp * Phi.diff
    grad$beta1[, j] = exp * (phi.diff / sqrt(h) + lcd$coef[2, j] * Phi.diff)
    grad$beta2[, j] = exp * (dphi.diff / h + 2 * lcd$coef[2, j] * phi.diff / sqrt(h) + lcd$coef[2, j]^2 * Phi.diff)
  }
  return(list(beta0 = apply(grad$beta0, 1, sum), beta1 = apply(grad$beta1, 1, sum), beta2 = apply(grad$beta2, 1, sum)))
}

line.slcd = function (lcd0, lcd1, x, w, h, ll0, tol=1e-10) {
  llt = function(alpha) {
    m = new.lcd((1 - alpha) * lcd0$alpha + alpha * lcd1$alpha, lcd1$theta,
                (1 - alpha) * lcd0$pi + alpha * lcd1$pi, lcd0$lower, lcd0$upper)
    eqx = Eqx(m, x, h)
    smo.f = eqx$x0[, 1] / m$C
    ll = w %*% (log(smo.f))
    list(lcd=m, smo.f=smo.f, ll=ll)
  }
  knots = c(lcd0$lower, lcd0$theta)
  grad = gradient.delta(lcd0, x, w, h, knots)
  delta = sum(grad * (c(lcd1$alpha, lcd1$pi) - c(lcd0$alpha, lcd0$pi))) * .333
  convergence = 0
  alpha = 1
  repeat{
    new = llt(alpha)
    if (new$ll >= ll0 + alpha * delta) break 
    if (alpha < tol) {
      convergence=1
      eqx = Eqx(lcd0, x, h)
      smo.f = eqx$x0[, 1] / lcd0$C
      new=list(lcd=lcd0, smo.f=smo.f, ll=ll0)
      break
    }
    alpha = alpha * .5
  }
  list(lcd=new$lcd, smo.f=new$smo.f, ll=new$ll, alpha=alpha, convergence=convergence)
}
