msvm_compact = function(K, y, lambda, epsilon = 1e-6, eig_tol_D = 0, epsilon_D = 1e-6)
{
  out = list()

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)
  n = length(y_int)
  qp_dim = n * n_class

  trans_Y = class_code(y_int, n_class)

  # Optimize alpha by solve.QP:
  # min (-d^Tb + 1/2 b^TDb)
  # subject to A^Tb <= b_0
  # Following steps (1) - (6)
  # (1) preliminary quantities

  Jk = matrix(1, nrow = n_class, ncol = n_class)
  Ik = diag(1, n_class)

  # vectorize y matrix
  y_vec = as.vector(trans_Y)

  # index for non-trivial alphas
  nonzeroIndex = y_vec != 1

  # (2) compute D
  D = (Ik - Jk / n_class) %x% K

  # subset the columns and rows for non-trivial alpha's
  Reduced_D = D[nonzeroIndex, nonzeroIndex]
  Reduced_D = fixit(Reduced_D, epsilon = eig_tol_D)
  max_D = max(abs(Reduced_D))
  diag(Reduced_D) = diag(Reduced_D) + max_D * epsilon_D

  # (3) compute g
  g = -n * lambda * y_vec

  # subset the components with non-trivial alpha's
  Reduced_g = g[nonzeroIndex]
  n_nonzeroIndex = length(Reduced_g)

  # equality constraint matrix
  R1 = ((Ik - Jk / n_class) %x% matrix(rep(1, n), nrow = 1))

  # eliminate one redundant equality constraint
  R1 = matrix(R1[1:(n_class - 1), ], nrow = n_class - 1, ncol = ncol(R1))

  # choose components with non-trivial alpha's
  Reduced_R1 = matrix(R1[, nonzeroIndex], nrow = nrow(R1), ncol = n_nonzeroIndex)

  # inequality constraint matrix
  R2 = diag(rep(1, n * (n_class - 1)))
  R2 = rbind(R2, -R2)

  # R consists of equality and inequality constraints
  R = t(rbind(Reduced_R1, R2))

  # (5) compute (b_0, b) = r
  # Right hand side of equality constraints
  r1 = rep(0, nrow(Reduced_R1))

  # Right hand side of inequality constraints
  r2 = c(rep(0, nrow(R2) / 2), rep(-1, nrow(R2) / 2))

  # R consists of right hand sides of equality and inequality constraints
  r = c(r1, r2)

  # (6) Find solution by solve.QP.compact
  nonzero = find_nonzero(R)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind
  dual = solve.QP.compact(Reduced_D, Reduced_g, Amat, Aind, r, meq = nrow(Reduced_R1))

  # Place the dual solution into the non-trivial alpha positions
  alpha = rep(0, qp_dim)
  alpha[nonzeroIndex] = dual$solution

  # Make alpha zero if they are too small
  alpha[alpha < 0] = 0
  # alpha[alpha > 1] = 1

  # Reshape alpha into a n by n_class matrix
  alpha = matrix(alpha, nrow = n)

  # Compute cmat = matrix of estimated coefficients
  cmat = -((alpha - matrix(rep(rowMeans(alpha), n_class), ncol = n_class)) / (n * lambda))

  # find b vector
  Kcmat = K %*% cmat

  # LP starts to find b vector
  # reformulate LP w/o equality constraint and redudancy
  # objective function with (b_j)_+,-(b_j)_, j=1,...,(k-1) and \xi_ij

  a = c(rep(0, 2 * (n_class - 1)), rep(1, n * (n_class - 1)))

  # inequality conditions
  B1 = -diag(1, n_class - 1) %x% rep(1, n)
  B2 = matrix(1, n, n_class - 1)
  A = cbind(B1, -B1)
  A = rbind(A, cbind(B2,-B2))
  A = A[nonzeroIndex, ] # reduced.A
  A = cbind(A, diag(1, n * (n_class - 1)))
  b = matrix(Kcmat - trans_Y, ncol = 1)
  b = b[nonzeroIndex] # reduced.b

  # constraint directions
  const.dir = matrix(rep(">=", nrow(A)))

  bpos = lp("min", objective.in = a, const.mat = A, const.dir = const.dir, const.rhs = b)$solution[1:(2 * (n_class - 1))]
  bvec = cbind(diag(1, n_class - 1), -diag(1, n_class - 1)) %*% matrix(bpos, ncol = 1)
  bvec = c(bvec, -sum(bvec))

  # Compute the fitted values
  fit = (matrix(rep(bvec, n), ncol = n_class, byrow = T) + Kcmat)
  fit_class = apply(fit, 1, which.max)

  # Return the output
  out$alpha = alpha
  out$cmat = cmat
  out$c0vec = bvec
  out$fit = fit
  out$fit_class = fit_class
  out$n = n
  out$n_class = n_class
  out$levels = levs
  return(out)
}

predict.msvm_compact = function(object, newK = NULL)
{
  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels
  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}

msvm = function(x = NULL, y, lambda, kernel, kparam, scale = FALSE, epsilon = 1e-6, eig_tol_D = 0, epsilon_D = 1e-8)
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

  K = kernelMatrix(x, x, kernel = kernel, kparam = kparam)
  solutions = msvm_compact(K = K, y = y, lambda = lambda, epsilon = epsilon, eig_tol_D = eig_tol_D, epsilon_D = epsilon_D)

  out$x = x
  out$y = y
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
  class(out) = "msvm"
  return(out)
}


predict.msvm = function(object, newx = NULL, newK = NULL)
{

  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }

  if (is.null(newK)) {
    newK = kernelMatrix(newx, object$x, kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels

  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}

cv.msvm = function(x, y, valid_x = NULL, valid_y = NULL, nfolds = 5, lambda_seq = 2^{seq(-10, 10, length.out = 100)},
                   kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
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
                            msvm_fit = msvm(x = x, y = y, lambda = params$lambda[j], kernel = kernel, kparam = params$kparam[j], scale = scale, ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.msvm(msvm_fit, newx = valid_x)$class

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
    opt_model = msvm(x = x, y = y, lambda = opt_param$lambda, kernel = kernel, kparam = opt_param$kparam, scale = scale, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "msvm"
  return(out)
}







































































