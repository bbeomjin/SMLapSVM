ramsvm_compact = function(K, y, gamma = 0.5, lambda, epsilon = 1e-6, eig_tol_D = 0, epsilon_D = 1e-6)
{
  out = list()

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(y_int)
  n = length(y_int)
  qp_dim = n * n_class

  code_mat = code_ramsvm(y_int)
  yyi = code_mat$yyi
  W = code_mat$W
  y_index = code_mat$y_index
  Hmatj = code_mat$Hmatj
  Lmatj = code_mat$Lmatj

  D = 0
  Amat = matrix(0, n * n_class, n_class - 1)
  for (j in 1:(n_class - 1)) {
    D = D + Hmatj[[j]] %*% K %*% t(Hmatj[[j]])
    Amat[, j] = -Lmatj[[j]]
  }

  D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(D))
  diag(D) = diag(D) + max_D * epsilon_D

  g_temp = matrix(-1, n, n_class)
  g_temp[y_index] = 1 - n_class
  g = as.vector(g_temp)

  dvec = -g * n * lambda

  Amat = cbind(Amat, diag(-1, n * n_class), diag(1, n * n_class))

  bvec_temp = matrix(gamma - 1, nrow = n, ncol = n_class)
  bvec_temp[y_index] = -gamma
  if (gamma == 0 | gamma == 1) {
    bvec_temp = bvec_temp - epsilon
  }
  # bvec = c(rep(0, qp_dim + n_class), as.vector(bvec_temp))
  bvec = c(rep(0, n_class - 1), as.vector(bvec_temp), rep(0, n * n_class))

  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, Amat, bvec, meq = n_class - 1)
  alpha = dual$solution
  alpha[alpha < 0] = 0

  alpha_mat = matrix(alpha, nrow = n, ncol = n_class)

  cmat = matrix(0, n, n_class - 1)
  for (j in 1:(n_class - 1)) {
    cmat[, j] = t(Hmatj[[j]]) %*% alpha / (n * lambda)
  }

  Kcmat = (K %*% cmat) %*% W

  alp_temp = matrix(1 - gamma, nrow = n, ncol = n_class)
  alp_temp[y_index] = gamma

  alp = c(as.vector(alp_temp), rep(0, 2 * (n_class - 1)))

  Alp1 = diag(qp_dim)
  Alp2 = matrix(0, nrow = qp_dim, ncol = 2 * (n_class - 1))

  for (i in 1:(n_class - 1)) {
    Alp2[, (2 * i - 1)] = Lmatj[[i]]
    Alp2[, (2 * i)] = -Lmatj[[i]]
  }

  Alp = cbind(Alp1, Alp2)

  blp_temp = Kcmat + 1
  blp_temp[y_index] = (n_class - 1) - Kcmat[y_index]
  blp = as.vector(blp_temp)

  const_dir = rep(">=", qp_dim)
  # const_dir[1] = "="
  cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir,const.rhs = blp)$solution[(qp_dim + 1):(qp_dim + 2 * (n_class - 1))]
  c0vec = rep(0, n_class - 1)
  for(j in 1:(n_class - 1)) {
    c0vec[j] = cposneg[(2 * j - 1)] - cposneg[(2 * j)]
  }

  W_c0vec = drop(t(c0vec) %*% W)
  # compute the fitted values
  fit = (matrix(W_c0vec, nrow = n, ncol = n_class, byrow = T) + Kcmat)
  fit_class = apply(fit, 1, which.max)
  # table(y, fit_class)

  # Return the output
  out$alpha = alpha_mat
  out$cmat = cmat
  out$c0vec = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n = n
  out$n_class = n_class
  out$levels = levs
  return(out)
}

predict.ramsvm_compact = function(object, newK = NULL) {

  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels
  n_class = object$n_class

  W = XI_gen(n_class)

  W_c0 = drop(t(c0vec) %*% W)

  pred_y = matrix(W_c0, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  # return(list(class = pred_y, inner_prod = inner_prod))
  return(list(class = pred_class, pred_value = pred_y))
}


ramsvm = function(x = NULL, y, gamma = 0.5, lambda, kernel, kparam, scale = FALSE, epsilon = 1e-6, eig_tol_D = 0, epsilon_D = 1e-8)
{
  out = list()
  n = NROW(x)
  p = ncol(x)

  center = rep(0, p)
  scaled = rep(1, p)

  if (scale) {
    x = scale(x)
    center = attr(x, "scaled:center")
    scaled = attr(x, "scaled:scale")
  }

  K = kernelMat(x, x, kernel = kernel, kparam = kparam)
  solutions = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda, epsilon = epsilon, eig_tol_D = eig_tol_D, epsilon_D = epsilon_D)

  out$x = x
  out$y = y
  out$gamma = gamma
  out$n_class = solutions$n_class
  out$levels = solutions$levels
  out$lambda = lambda
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$c0vec = solutions$c0vec
  out$alpha = solutions$alpha
  out$epsilon = epsilon
  out$eig_tol_D = eig_tol_D
  out$epsilon_D = epsilon_D
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "ramsvm"
  return(out)
}


predict.ramsvm = function(object, newx = NULL, newK = NULL, ...) {

  if (is.null(newK)) {
    newK = kernelMat(newx, object$x, kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels
  n_class = object$n_class

  W = XI_gen(n_class)

  W_c0 = drop(t(c0vec) %*% W)

  pred_y = matrix(W_c0, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_class = levs[apply(fit, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}




cv.ramsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5, lambda_seq = 2^{seq(-10, 10, length.out = 100)},
                      kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                      scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  # if (scale) {
  #   x = scale(x)
  #   if (!is.null(valid_x)) {
  #     means = attr(x, "scaled:center")
  #     stds = attr(x, "scaled:scale")
  #     valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE)) / matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
  #   }
  # }

  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)

  # The number of classes
  n_class = length(unique(y))

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  kparam = sort(kparam, decreasing = TRUE)

  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, kparam = kparam)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            msvm_fit = ramsvm(x = x, y = y, gamma = gamma, lambda = params$lambda[j], kernel = kernel, kparam = params$kparam[j], scale = scale, ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.ramsvm(msvm_fit, newx = valid_x)$class

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val) / length(valid_y)
                              err = 1 - acc
                            } else {
                              # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                          } else {
                            msvm_fit = NULL
                            err = Inf
                          }

                          return(list(error = err, fit_model = msvm_fit))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }

  out$opt_param = c(lambda = opt_param$lambda, kparam = opt_param$kparam)
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$scale = scale
  if (optModel) {
    opt_model = ramsvm(x = x, y = y, gamma = gamma, lambda = opt_param$lambda, kernel = kernel, kparam = opt_param$kparam, scale = scale, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "ramsvm"
  return(out)
}








