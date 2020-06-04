# Smoothed Log-concave Regression
slcr = function (y, z, beta=NULL, lcd, h, tol=1e-10, maxit=100) {
  z = as.matrix(z)
  # The number of observations and variables
  n = length(y)
  p = ifelse(is.na(nrow(z)), 1, ncol(z))
  # Initial regression parameter
  if (is.null(beta)) beta = lm(y ~ 0+z)$coefficients
  # Initial log-concave error density
  if (p == 1) e = y - z*beta else e = y - z%*%beta
  if (missing(lcd)) lcd = cnmlcd(e)$lcd
  if (missing(h)) {
    nk = length(lcd$theta) + 2
    h = var(x) + mean(x)^2 - lcd$cpk[3, nk - 1]
  }
  res = cnmslcd(x=e, lcd=lcd, h=h)
  lcd = res$lcd
  smo.f = res$smo.f
  ll = res$ll
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    beta.old = beta
    lcd.old = lcd
    smo.f.old = smo.f
    ll.old = ll
    # Step 1: Update regression parameter estimate
    ll.temp = -Inf
    for (i in 1:maxit) {
      if (ll <= ll.temp + tol) break
      beta.temp = beta
      ll.temp = ll
      grad = gradient.beta(lcd, e, h=h)
      W1 = -grad$beta1 / grad$beta0
      W2 = diag((grad$beta0 * grad$beta2 - grad$beta1^2) / grad$beta0^2)
      beta = drop(solve(t(z)%*%W2%*%z)%*%t(z)%*%(W2%*%z%*%beta - W1))
      if (p == 1) e = y - z*beta else e = y - z%*%beta
      ll = logLik.slcd(lcd, e, h)
    }
    # Step 2: Update log-concave error density estimate
    res = cnmslcd(x=e, h=h)
    lcd = res$lcd
    smo.f = res$smo.f
    ll = res$ll
  }
  list(beta=beta, lcd=res$lcd, ll=ll, num.iterations=i, convergence=convergence)
}

# Doubly Smoothed Log-concave Regression
dslcr = function (y, z, beta=NULL, h, tol=1e-10, maxit=100) {
  z = as.matrix(z)
  y = as.matrix(y)
  z = rbind(z, z, z, z, z)
  qn = qnorm((1:5)/6, mean=0, sd=sqrt(h))
  qn = qn / sd(qn)
  y = rbind(y + qn[1], y + qn[2], y + qn[3], y + qn[4], y + qn[5])
  # The number of observations and variables
  n = length(y)
  p = ifelse(is.na(nrow(z)), 1, ncol(z))
  # Initial regression parameter
  if (is.null(beta)) beta = lm(y ~ 0+z)$coefficients
  # Initial log-concave error density
  if (p == 1) e = y - z*beta else e = y - z%*%beta
  if (missing(lcd)) lcd = cnmlcd(e)$lcd
  if (missing(h)) {
    nk = length(lcd$theta) + 2
    h = var(x) + mean(x)^2 - lcd$cpk[3, nk - 1]
  }
  res = cnmslcd(x=e, lcd=lcd, h=h)
  lcd = res$lcd
  smo.f = res$smo.f
  ll = res$ll
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    beta.old = beta
    lcd.old = lcd
    smo.f.old = smo.f
    ll.old = ll
    # Step 1: Update regression parameter estimate
    ll.temp = -Inf
    for (i in 1:maxit) {
      if (ll <= ll.temp + tol) break
      beta.temp = beta
      ll.temp = ll
      grad = gradient.beta(lcd, e, h=h)
      W1 = -grad$beta1 / grad$beta0
      W2 = diag((grad$beta0 * grad$beta2 - grad$beta1^2) / grad$beta0^2)
      beta = drop(solve(t(z)%*%W2%*%z)%*%t(z)%*%(W2%*%z%*%beta - W1))
      if (p == 1) e = y - z*beta else e = y - z%*%beta
      ll = logLik.slcd(lcd, e, h)
    }
    # Step 2: Update log-concave error density estimate
    res = cnmslcd(x=e, h=h)
    lcd = res$lcd
    smo.f = res$smo.f
    ll = res$ll
  }
  list(beta=beta, lcd=res$lcd, ll=ll, num.iterations=i, convergence=convergence)
}
