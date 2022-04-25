generateMultiorange = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1)
{
  set.seed(seed)
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 1.5)
    sx = sum(x^2)
    if (sx <= 0.5) {
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }
    else if (1.5 < sx & sx <= 2.5) {
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
    else if (4.5 < sx & sx <= 6.5) {
      y[k] = 3
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = 1.5), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}

generateMultiorange2 = function(n, p = 2, seed = 1, with_noise = TRUE, noise_p = 1)
{
  set.seed(seed)
  X = matrix(nrow = n, ncol = p)
  y = numeric(n)
  k = 1
  while (k <= n) {
    x = rnorm(p, sd = 2)
    sx = sum(x^2)
    if (sx <= 0.8) {
      y[k] = 1
      X[k, ] = x
      k = k + 1
    }
    else if (2.5 < sx & sx <= 4) {
      y[k] = 2
      X[k, ] = x
      k = k + 1
    }
    else if (7 < sx & sx <= 9.5) {
      y[k] = 3
      X[k, ] = x
      k = k + 1
    }
  }
  if (with_noise) {
    noise_dat = matrix(rnorm(n * noise_p, sd = 2), n, noise_p)
    X = cbind(X, noise_dat)
  }
  return(list(x = X, y = y))
}


generateMultiMoon = function(each_n = 100, sigma = 1, noise_p = 4, noise_sd = 3, seed = NULL)
{
  set.seed(seed)
  x = runif(each_n, 0, pi)
  c1 = cbind(5 * cos(x) - 3.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, pi, 2 * pi)
  c2 = cbind(5 * cos(x) + 3.5 + rnorm(each_n) * sigma, 10 * sin(x) +
               2.5 + rnorm(each_n) * sigma)
  x = runif(each_n, 0, pi)
  c3 = cbind(5 * cos(x) + 10.5 + rnorm(each_n) * sigma, 10 * sin(x) -
               2.5 + rnorm(each_n) * sigma)
  X = rbind(c1, c2, c3)
  noise_X = matrix(rnorm(3 * each_n * noise_p, 0, noise_sd), nrow = 3 * each_n, ncol = noise_p)
  X = cbind(X, noise_X)
  y = rep(c(1, 2, 3), each = each_n)
  return(list(x = X, y = y))
}


sim_gen = function(n, p, seed = NULL, type = c("bayes", "poly", "cosso", "neuralnet"))
{
  call = match.call()
  type = match.arg(type)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (type == "cosso") {
    X = matrix(rnorm(n * p, 0, sd = 2), n, p)

    g1 = function(x) {x}
    g2 = function(x) {(2 * x - 1)^2}


    lr1 = 1.1 * g1(X[, 1]) + 1.6 * g2(X[, 2]) - 8.2 #+ 1 * g3(X[, 3]) - 6 * g4(X[, 4])
    lr2 = 3.3 * g1(X[, 1]) + 1.5 * g2(X[, 2]) - 7.0 #- 3.5 * g3(X[, 3]) - 4.5 * g4(X[, 4])

    const = (1 + exp(lr1) + exp(lr2))
    prob1 = exp(lr1) / const
    prob2 = exp(lr2) / const
    prob3 = 1 / const
    probs = cbind(prob1, prob2, prob3)

    # y = apply(probs, 1, function(prob) {
    #   sample(1:3, 1, TRUE, prob)})

    y = apply(probs, 1, which.max)

    out = list()
    out$x = X
    out$y = y
    out$true = rep(c(1, 0), c(2, p - 2))
  }

  if (type == "bayes") {
    dat = mlbench.2dnormals(n = n, cl = 3, sd = 1)
    X_true = dat$x
    r = ncol(X_true)
    y = dat$classes
    X_noise = matrix(rnorm(n * (p - r), sd = 1), nrow = n, ncol = p - r)
    X = cbind(X_true, X_noise)
    out = list()
    out$x = X
    out$y = y
    out$true = rep(c(1, 0), c(r, p - r))
  }

  if (type == "poly") {
    r = 2
    X_tmp = matrix(rnorm(n * r, 0, 0.8), n, r)
    x1 = X_tmp[, 1]; x2 = X_tmp[, 2]
    c = 1
    X_kern = data.matrix(data.frame(x1^3, x2^3, sqrt(3) * x1^2 * x2, sqrt(3) * x1 * x2^2,
                                    sqrt(3 * c) * x1^2, sqrt(3 * c) * x2^2, sqrt(6 * c) * x1 * x2,
                                    sqrt(3) * c * x1, sqrt(3) * c * x2, sqrt(c^3)))
    beta1 = c(3, -2, 3.5, -3.5, 5.5, -4.5, 3, 0, -2, 0)
    beta2 = c(3, -3, 2.5, -2.5, 1, -1, 0, -2, -2, 0)
    beta3 = beta1 - beta2
    lr = X_kern %*% cbind(beta1, beta2, beta3)

    probs = exp(lr - as.vector(HTLR:::log_sum_exp(lr)))
    # y = apply(probs, 1, function(prob) {
    #   sample(1:3, 1, TRUE, prob)})
    y = apply(probs, 1, which.max)
    X_noise = matrix(rnorm(n * (p - r)), n, (p - r))
    X = cbind(X_tmp, X_noise)
    out = list()
    out$x = X
    out$y = y
    out$true = rep(c(1, 0), c(r, p -r))
  }

  if (type == "neuralnet") {
    r = 2
    X_tmp = matrix(rnorm(n * r, 0, 2.0), n, r)
    X_tmp = cbind(X_tmp, X_tmp[, 1] * X_tmp[, 2])
    sigmoid = function(x) {
      return(1 / (1 + exp(-x)))
    }

    LeLU = function(x) {
      return(pmax(0, x))
    }

    node = c(4, 3, 2)

    beta_mat1 = matrix(c(3.5, 1.3, -1.5, 2.3,
                         0.5, 3.2, 5.5, -3.0,
                         -3.0, -4.2, -3.0, -1.0),
                       nrow = 3, byrow = TRUE)

    node1_mat = drop(X_tmp %*% beta_mat1) + 1
    layer1 = matrix(sigmoid(node1_mat), nrow = n, ncol = node[1])

    beta_mat2 = matrix(c(2.8, -2.5, -1.5,
                         1.2, -2.0, -3.0,
                         -2.7, 1.0, 2.0,
                         -3.0, 2.0, 1.0),
                       nrow = node[1], ncol = node[2], byrow = TRUE)
    node2_mat = drop(layer1 %*% beta_mat2) - 1
    layer2 = matrix(sigmoid(node2_mat), nrow = n, ncol = node[2])

    output_beta1 = c(3.1, -2.9, 0.5)
    output_beta2 = c(-1.7, 4.1, 1.6)
    output_beta3 = c(2.5, -2, 4.1)

    prob1 = sigmoid(drop(layer2 %*% output_beta1))
    prob2 = sigmoid(drop(layer2 %*% output_beta2))
    prob3 = sigmoid(drop(layer2 %*% output_beta3))
    probs = cbind(prob1, prob2, prob3)

    y = apply(probs, 1, which.max)

    noise = matrix(rnorm(n * (p - r), 0, 2.0), n, (p - r))
    X = cbind(X_tmp[, -3], noise)
    true_vec = c(rep(1, r), rep(0, p - r))
    out = list()
    out$x = X
    out$y = y
    out$true = true_vec
    out$probs = probs
  }

  return(out)
}


