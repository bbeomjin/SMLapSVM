# dyn.load("../src/alpha_update.dll")
ramlapsvm_core = function(K, L, y, gamma = 0.5, lambda, lambda_I, weight = NULL, epsilon = 1e-4 * length(y) * length(unique(y)), maxiter = 300)
{

  out = list()
  n_class = length(unique(y))
  n_l = length(y)
  n = nrow(K)
  n_u = n - n_l

  inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  Q = (n_l * lambda) * (K %*% inv_LK)[1:n_l, 1:n_l] + 1

  warm = matrix(data = 0.0, nrow = n_l, ncol = n_class)
  if (is.null(weight)) {weight = numeric(n_l) + 1.0}

  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  my = t(XI_gen(k = n_class))

  yyi = Y_matrix_gen(k = n_class, nobs = n_l, y = y)

  alpha_ij = warm
  alpha_yi = numeric(n_l)


  erci = as.double(-rep(1, ncol(Q)) / n_l / lambda)

  aa = .C("alpha_update",
            as.vector(alpha_ij),
            as.vector(alpha_yi),
            as.vector(my),
            as.vector(yyi),
            as.vector(Q),
            as.double(lambda),
            as.vector(weight),
            as.integer(n_l),
            as.double(n_l),
            as.integer(n_class),
            as.double(n_class),
            as.vector(erci),
            as.double(gamma),
            as.vector(as.integer(y)),
            as.double(epsilon),
            outalpha_ij = as.vector(numeric(n_l * n_class)),
            maxiter = as.integer(maxiter), PACKAGE = "SMLapSVM")

  warm = matrix(data = aa$outalpha_ij, nrow = n_l, ncol = n_class)

  beta = beta_kernel(y = y,
                     k = n_class,
                     my = my,
                     warm = warm,
                     lambda = lambda,
                     inv_LK = inv_LK)
  # drop(crossprod(beta[[1]][, 1], X))

  # tt = beta_linear(x = x, y = y_train, k = k, my = my, warm = warm, lambda = templambda)

  # beta0 = matrix(find_theta2(y, K, gamma = gamma, cmat = beta$beta, lambda = lambda), nrow = 1)

  betaout = beta$beta
  beta0out = beta$beta0
  # beta0out[[count]] = beta0

  out$y = y
  out$n_class = n_class
  out$gamma = gamma
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$beta = betaout
  out$beta0 = beta0out
  out$epsilon = epsilon
  out$warm = warm
  class(out) = "ramlapsvm_core"
  return(out)
}

predict.ramlapsvm_core = function(object, newK = NULL) {

  beta = object$beta
  beta0 = object$beta0

  temp_pred_y = predict_kernel(K_test = newK,
                               beta = beta,
                               beta0 = beta0,
                               k = object$n_class)
  inner_prod = temp_pred_y$inner_prod
  temp_pred_y = temp_pred_y$class
  pred_y = temp_pred_y

  return(list(class = pred_y, inner_prod = inner_prod))
}


