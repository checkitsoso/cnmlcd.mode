cnmslcd = function (x, lcd, h, maxit=100, tol=1e-06) {
  require(lsei)
  lambda = 1e-15
  xw = x.weight(x)
  x = xw$x
  w = xw$w
  n = sum(w)
  nx = length(x)
  lower = x[1]
  upper = x[nx]
  if (missing(lcd)) lcd = cnmlcd(x)$lcd
  if (missing(h)) {
    nk = length(lcd$theta) + 2
    h = var(x) + mean(x)^2 - lcd$cpk[3, nk - 1]
  }
  nc = dim(lcd$coef)[2]
  eqx = Eqx(lcd, x, h)
  smo.f = eqx$x0[, 1] / lcd$C
  ll = w %*% log(smo.f)
  delta.cand = sort(c(lcd$theta, seq(lower, upper, length=200)))
  convergence = 1
  ll.old = -Inf
  for (i in 1:maxit) {
    if (ll <= ll.old + tol) {
      convergence = 0
      break
    }
    lcd.old = lcd
    ll.old = ll
    gd = maxima.delta(lcd, x, w, h, delta.cand)
    # plot(delta.cand, gd$gds, type="l")
    # abline(h=0, lty=2, col="red")
    if (length(gd$delta) != 0) {
      nsp = gd$delta
      nsl = length(nsp)
      if (nsl >= 1) lcd = new.lcd(lcd$alpha, c(lcd$theta, nsp), c(lcd$pi, double(nsl)), lcd$lower, lcd$upper)
    }
    lcd1 = lcd
    knots = c(lcd1$lower, lcd1$theta, lcd1$upper)
    for (j in 1:maxit) {
      gdk = gradient.delta(lcd1, x, w, h, knots)
      if (all(gdk < tol) | any(is.nan(gdk))) break
      gt = gradient.theta(lcd1, x, w, h)
      e = eigen(gt$H)
      v2 = sqrt(e$values[e$values >= e$values[1] * lambda])
      kr = length(v2)
      R = t(e$vectors[, 1:kr, drop = FALSE]) * v2
      p = gt$grad/n + drop(c(-lcd1$alpha, lcd1$pi) %*% gt$H)
      b = drop(crossprod(e$vectors[, 1:kr, drop = FALSE], p))/v2
      r1 = lsei::pnnls(R, b, 1)
      lcd1 = new.lcd(alpha = -r1$x[1], theta = lcd1$theta, pi = r1$x[-1], lower = lower, upper = upper)
    }
    r = line.slcd(lcd, lcd1, x, w, h, ll.old)
    lcd = r$lcd
    ll = r$ll
    if (any(lcd$pi == 0)) lcd = simplify.lcd(lcd)
  }
  list(lcd = lcd, smo.f = r$smo.f, ll = ll, num.iterations = i, max.gradient = gd$gmax, convergence = convergence)
}

