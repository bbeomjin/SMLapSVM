rmlapsvm_compact = function(K, L, y, gamma = 0.5, lambda, lambda_I, epsilon = 1e-6, eig_tol_D = 1e-15, eig_tol_I = 2e-15)
{
  out = list()
  # The labeled sample size, unlabeled sample size, the number of classes and dimension of QP problem
  n_class = length(unique(y))
  n = nrow(K)
  n_l = length(y)
  n_u = n - n_l
  qp_dim = n_l * n_class

  code_mat = code(y)
  In = code_mat$In
  vmatj = code_mat$vmatj
  umatj = code_mat$umatj
  Hmatj = code_mat$Hmatj
  y_index = code_mat$y_index

  J = cbind(diag(n_l), matrix(0, n_l, n_u))
  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  LK = fixit(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K), eig_tol_I)
  inv_LK = chol2inv(chol(LK))
  Q = J %*% K %*% inv_LK %*% t(J)
  # Q = J %*% t(inv_LK) %*% K %*% t(J)

  # Compute Q = K x inv_LK
  D = matrix(0, qp_dim, qp_dim)
  Amat = matrix(0, (2 * qp_dim + n_class), qp_dim)

  for (k in 1:n_class) {
    D = D + t(Hmatj[[k]]) %*% Q %*% Hmatj[[k]]
    Amat[k, ] = rep(1, n_l) %*% Hmatj[[k]]
  }

  # max_D = max(abs(D))
  # D = D / max_D
  D = fixit(D, epsilon = eig_tol_D)
  # diag(D) = diag(D) + epsilon_D

  g_temp = matrix(-1, n_l, n_class)
  g_temp[y_index] = -n_class + 1
  g = as.vector(g_temp)

  # g = rep(-1, qp_dim)
  # for(j in 1:n_class) {
  #   for(i in 1:n_l) {
  #     if (y[i] == j) {
  #       g[(j - 1) * n_l + i] = -(n_class - 1)
  #     }
  #   }
  # }

  dvec = -g
  # dvec = -g / max_D

  diag(Amat[(n_class + 1):(n_class + qp_dim), ]) = 1
  diag(Amat[(n_class + qp_dim + 1):(n_class + 2 * qp_dim), ]) = -1

  # (3) compute Ama

  # (4) compute bvec
  # bvec = rep(0, (2 * qp_dim + n_class))

  bvec_temp = matrix(gamma - 1, nrow = n_l, ncol = n_class)
  bvec_temp[y_index] = -gamma
  if (gamma == 0 | gamma == 1) {
    bvec_temp = bvec_temp - epsilon
  }
  bvec = c(rep(0, qp_dim + n_class), as.vector(bvec_temp))

  # for (j in 1:n_class) {
  #   for (i in 1:n_l) {
  #     flag = 0
  #     if (y[i] == j) {
  #       flag = 1
  #     }
  #     bvec[n_class + qp_dim + (j - 1) * n_l + i] = -(gamma * flag + (1 - gamma) * (1 - flag))
  #     # correction to avoid redundant constraints when gamma = 0 or 1
  #     if ((flag == 1 & gamma == 0) | (flag == 0 & gamma == 1)) {
  #       bvec[n_class + qp_dim + (j - 1) * n_l + i] = bvec[n_class + qp_dim + (j - 1) * n_l + i] - epsilon_D
  #     }
  #   }
  # }

  # remove one redudant constraint
  Amat1 = Amat[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class)), ]
  bvec1 = bvec[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class))]

  # (5) find solution by solve.QP

  nonzero = find_nonzero(t(Amat1))
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec1, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, t(Amat1), bvec1, meq = (n_class - 1))

  alpha = dual$solution
  alpha[alpha < 0] = 0

  alpha_mat = matrix(alpha, nrow = n_l, ncol = n_class)
  alpha_mat[y_index][alpha_mat[y_index] > gamma] = gamma

  for (j in 1:n_class) {
    alpha_mat[y != j, j][alpha_mat[y != j, j] > (1 - gamma)] = (1 - gamma)
  }

  # for (j in 1:n_class) {
  #   for (i in 1:n_l) {
  #     if (y[i] == j & (alpha[(j - 1) * n_l + i] > gamma)) {
  #       alpha[(j - 1) * n_l + i] = gamma
  #     }
  #     if (y[i] != j & (alpha[(j - 1) * n_l + i] > (1 - gamma))) {
  #       alpha[(j - 1) * n_l + i] = (1 - gamma)
  #     }
  #   }
  # }

  cmat_temp = matrix(0, n_l, n_class)
  for (k in 1:n_class) {
    cmat_temp[, k] = Hmatj[[k]] %*% alpha
  }
  cmat = inv_LK %*% t(J) %*% cmat_temp

  # find b vector using LP
  Kcmat = J %*% K %*% cmat

  alp_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  alp_temp[y_index] = gamma

  alp = c(as.vector(alp_temp), rep(0, 2 * n_class))

  # alp = rep((1 - gamma), (qp_dim + 2 * n_class))
  # for (j in 1:n_class) {
  #   for (i in 1:n_l) {
  #     if (y[i] == j) {
  #       alp[n_l * (j - 1) + i] = gamma
  #     }
  #   }
  # }
  # alp[(qp_dim + 1):(qp_dim + 2 * n_class)] = 0

  # constraint matrix and vector

  Alp1 = c(rep(0, qp_dim), rep(c(1, -1), n_class))
  Alp2 = diag(qp_dim)

  Alp3 = matrix(0, nrow = qp_dim, ncol = 2 * n_class)

  Alp3_temp = matrix(-1, nrow = n_l, ncol = n_class)
  Alp3_temp[y_index] = 1

  for (i in 1:n_class) {
    Alp3[(n_l * (i - 1) + 1):(n_l * i), (2 * i - 1)] = Alp3_temp[, i]
    Alp3[(n_l * (i - 1) + 1):(n_l * i), (2 * i)] = -Alp3_temp[, i]
  }

  Alp = rbind(Alp1, cbind(Alp2, Alp3))

  blp_temp = Kcmat + 1
  blp_temp[y_index] = (k - 1) - Kcmat[y_index]
  blp = c(0, as.vector(blp_temp))

  # constraint directions
  const_dir = rep(">=", (qp_dim + 1))
  const_dir[1] = "="
  cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir,const.rhs = blp)$solution[(qp_dim + 1):(qp_dim + 2 * n_class)]
  c0vec = rep(0, n_class)
  for(j in 1:n_class) {
    c0vec[j] = cposneg[(2 * j - 1)] - cposneg[(2 * j)]
  }

  # compute the fitted values
  fit = (matrix(rep(c0vec, n_l), ncol = n_class, byrow = T) + Kcmat)
  fit_class = apply(fit, 1, which.max)

  # Return the output
  out$alpha = alpha_mat
  out$cmat = cmat
  out$c0vec = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n_l = n_l
  out$n_u = n_u
  out$n_class = n_class
  return(out)
}


rmlapsvm = function(x = NULL, y = NULL, ux = NULL, gamma = 0.5, lambda, lambda_I, kernel, kparam, scale = FALSE,
                    adjacency_k = 6, normalized = TRUE, weight = NULL, weightType = "Binary", epsilon = 1e-6, eig_tol_D = 1e-15, eig_tol_I = 2e-15)
{
  out = list()
  n_l = NROW(x)
  n_u = NROW(ux)
  rx = rbind(x, ux)
  p = ncol(x)
  center = rep(0, p)
  scaled = rep(1, p)

  if (scale) {
    rx = scale(rx)
    center = attr(rx, "scaled:center")
    scaled = attr(rx, "scaled:scale")
    x = (x - matrix(center, nrow = n_l, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_l, ncol = p, byrow = TRUE)
    ux = (ux - matrix(center, nrow = n_u, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_u, ncol = p, byrow = TRUE)
  }

  n = n_l + n_u
  n_class = max(y)

  K = kernelMat(rx, rx, kernel = kernel, kparam = kparam)

  # W = RSSL:::adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # d = rowSums(W)
  # L = diag(d) - W
  # if (normalized) {
  #   L = diag(1 / sqrt(d)) %*% L %*% diag(1 / sqrt(d))
  # }



  # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # graph = W
  graph = make_knn_graph_mat(rx, k = adjacency_k)
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)

  solutions = rmlapsvm_compact(K = K, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I,
                               epsilon = epsilon, eig_tol_D = eig_tol_D, eig_tol_I = eig_tol_I)

  out$x = x
  out$ux = ux
  out$y = y
  out$n_class = n_class
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$c0vec = solutions$c0vec
  out$alpha = solutions$alpha
  out$epsilon = epsilon
  out$eig_tol_D = eig_tol_D
  out$eig_tol_I = eig_tol_I
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "rmlapsvm"
  return(out)
}

predict.rmlapsvm_compact = function(object, newK = NULL)
{
  cmat = object$cmat
  c0vec = object$c0vec
  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = apply(pred_y, 1, which.max)
  return(list(class = pred_class, pred_value = pred_y))
}

predict.rmlapsvm = function(object, newx = NULL, newK = NULL)
{

  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }

  if (is.null(newK)) {
    newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = object$cmat
  c0vec = object$c0vec

  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = apply(pred_y, 1, which.max)
  return(list(class = pred_class, pred_value = pred_y))
}

Kfold_rmlapsvm = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                         gamma = 0.5, lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                         kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1), normalized = FALSE,
                         scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
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

  if (!is.numeric(lambda_seq)) {
    lambda_seq = as.numeric(lambda_seq)
  }

  if (!is.numeric(lambda_I_seq)) {
    lambda_I_seq = as.numeric(lambda_I_seq)
  }

  if (!is.numeric(kparam)) {
    kparam = as.numeric(kparam)
  }

  # The number of classes
  k = length(unique(y))

  # Combination of hyper-parameters
  # params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq[order(lambda_I_seq, decreasing = TRUE)], kparam = kparam)
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq, kparam = kparam)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            rmsvm_fit = rmlapsvm(x = x, y = y, ux = ux, gamma = gamma, lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                                 kernel = kernel, kparam = params$kparam[j], scale = scale, normalized = normalized, ...)
                          })

                          if (!inherits(error, "try-error")){
                            pred_val = predict.rmlapsvm(rmsvm_fit, newx = valid_x)$class

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val) / length(valid_y)
                              err = 1 - acc
                            } else {
                              # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                          } else {
                            rmsvm_fit = NULL
                            err = 1
                          }

                          return(list(error = err, fit_model = rmsvm_fit))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }

  out = list()
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$fold_models = lapply(model_list, "[[", opt_ind)
  out$fold_ind = fold_list
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$scale = scale
  if (optModel) {
    opt_model = rmlapsvm(x = x, y = y, ux = ux, gamma = gamma, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
                        kernel = kernel, kparam = opt_param$kparam, scale = scale, normalized = normalized, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "rmlapsvm"
  return(out)
}
