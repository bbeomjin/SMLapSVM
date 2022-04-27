rmlapsvm_compact = function(K, L, y, gamma = 0.5, lambda, lambda_I, epsilon = 1e-6,
                            eig_tol_I = 1e-13,
                            eig_tol_D = 0,
                            inv_tol = 1e-25, epsilon_D = 1e-8)
{
  out = list()
  # The labeled sample size, unlabeled sample size, the number of classes and dimension of QP problem
  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)

  n = nrow(K)
  n_l = length(y_int)
  n_u = n - n_l
  qp_dim = n_l * n_class

  code_mat = code_rmsvm(y_int)
  In = code_mat$In
  # vmatj = code_mat$vmatj
  # umatj = code_mat$umatj
  Hmatj = code_mat$Hmatj
  y_index = code_mat$y_index

  J = cbind(diag(n_l), matrix(0, n_l, n_u))
  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  # LK = fixit(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K), inv_tol)
  # inv_LK = chol2inv(chol(LK))
  # K = fixit(K, eig_tol_D)
  # inv_LK = inverse(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K), epsilon = inv_tol)
  LK = diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K)
  # LK = fixit(LK, epsilon = eig_tol_I)
  # max_LK = max(abs(LK))
  # inv_LK = solve(LK + diag(max_LK * epsilon_I, n), t(J))

  JK = J %*% K

  # inv_LK = solve(LK, t(J), tol = inv_tol)
  inv_LK = solve(LK, tol = inv_tol) %*% t(J)
  # inv_LK = solve(LK, t(J))
  # inv_LK = solve(LK + diag(max(abs(diag(LK))) * epsilon_I, n), tol = inv_tol)
  # inv_LK = inv_LK %*% t(J)
  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))

  # inv_LK = solve(LK / max_LK + diag(epsilon_I, n), t(J) / max_LK)

  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))
  # inv_LK = inverse(LK, epsilon = inv_tol)

  Q = JK %*% inv_LK
  # Q = (Q + t(Q)) / 2
  # Q = J %*% K %*% inv_LK %*% t(J)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # Q = J %*% t(inv_LK) %*% K %*% t(J)

  # Compute Q = K x inv_LK
  D = matrix(0, qp_dim, qp_dim)
  Amat = matrix(0, (2 * qp_dim + n_class), qp_dim)
  for (j in 1:n_class) {
    D = D + t(Hmatj[[j]]) %*% Q %*% Hmatj[[j]]
    Amat[j, ] = rep(1, n_l) %*% Hmatj[[j]]
  }
  D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(diag(D)))
  # D = D / max_D
  diag(D) = diag(D) + max_D * epsilon_D

  # diag(D) = diag(D) + epsilon_D

  g_temp = matrix(-1, n_l, n_class)
  g_temp[y_index] = -n_class + 1
  g = as.vector(g_temp)

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
  Amat = Amat[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class)), ]
  bvec = bvec[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class))]

  # (5) find solution by solve.QP

  nonzero = find_nonzero(t(Amat))
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, t(Amat1), bvec1, meq = (n_class - 1))

  alpha = dual$solution
  alpha[alpha < 0] = 0

  alpha_mat = matrix(alpha, nrow = n_l, ncol = n_class)
  # alpha_mat[y_index][alpha_mat[y_index] > gamma] = gamma

  # for (j in 1:n_class) {
  #   alpha_mat[y != j, j][alpha_mat[y != j, j] > (1 - gamma)] = (1 - gamma)
  # }
  #
  # alpha = as.vector(alpha_mat)
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

  cmat = matrix(0, n, n_class)
  for (k in 1:n_class) {
    cmat[, k] = inv_LK %*% Hmatj[[k]] %*% alpha
  }

  # find b vector using LP
  Kcmat = JK %*% cmat

  upper_mat = matrix(1 - gamma, nrow = nrow(alpha_mat), ncol = n_class)
  upper_mat[y_index] = gamma
  logic = (alpha_mat > 0) & (alpha_mat < upper_mat)
  if (all(colSums(logic) > 0)) {
    c0mat = -Kcmat - 1
    c0mat[y_index] = (n_class - 1) - Kcmat[y_index]
    c0vec = colMeans(c0mat)
  } else {
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
    blp_temp[y_index] = (n_class - 1) - Kcmat[y_index]
    blp = c(0, as.vector(blp_temp))

    # constraint directions
    const_dir = rep(">=", (qp_dim + 1))
    const_dir[1] = "="
    cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir,const.rhs = blp)$solution[(qp_dim + 1):(qp_dim + 2 * n_class)]
    c0vec = rep(0, n_class)
    for(j in 1:n_class) {
      c0vec[j] = cposneg[(2 * j - 1)] - cposneg[(2 * j)]
    }
  }

  # compute the fitted values
  fit = (matrix(rep(c0vec, n_l), ncol = n_class, byrow = T) + Kcmat)
  fit_class = levs[apply(fit, 1, which.max)]
  if (attr(levs, "type") == "factor") {fit_class = factor(fit_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {fit_class = as.numeric(fit_class)}
  if (attr(levs, "type") == "integer") {fit_class = as.integer(fit_class)}

  # Return the output
  out$alpha = alpha_mat
  out$cmat = cmat
  out$c0vec = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n_l = n_l
  out$n_u = n_u
  out$n_class = n_class
  out$levels = levs
  return(out)
}


predict.rmlapsvm_compact = function(object, newK = NULL)
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


rmlapsvm = function(x = NULL, y = NULL, ux = NULL, gamma = 0.5, lambda, lambda_I, kernel, kparam, scale = FALSE,
                    adjacency_k = 6, normalized = TRUE, weight = NULL, weightType = "Binary", epsilon = 1e-6,
                    eig_tol_I = 1e-13, eig_tol_D = 0, inv_tol = 1e-25, epsilon_D = 1e-8)
{
  out = list()
  n_l = NROW(x)
  n_u = NROW(ux)
  n = n_l + n_u
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


  K = kernelMatrix(rx, rx, kernel = kernel, kparam = kparam)

  # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # graph = W
  graph = make_knn_graph_mat(rx, k = adjacency_k)
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)

  solutions = rmlapsvm_compact(K = K, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, epsilon = epsilon,
                               eig_tol_D = eig_tol_D, inv_tol = inv_tol, epsilon_D = epsilon_D)

  out$x = x
  out$ux = ux
  out$y = y
  out$n_class = solutions$n_class
  out$levels = solutions$levels
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$c0vec = solutions$c0vec
  out$alpha = solutions$alpha
  out$epsilon = epsilon
  out$eig_tol_D = eig_tol_D
  out$inv_tol = inv_tol
  out$epsilon_D = epsilon_D
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "rmlapsvm"
  return(out)
}


predict.rmlapsvm = function(object, newx = NULL, newK = NULL)
{

  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }

  if (is.null(newK)) {
    newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
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

cv.rmlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                       lambda_seq = 2^{seq(-10, 15, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                       kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                       scale = FALSE, adjacency_k = 6, weightType = "Heatmap", normalized = TRUE,
                       criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
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
  lambda_I_seq = as.numeric(lambda_I_seq)
  kparam = as.numeric(kparam)

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  lambda_I_seq = sort(lambda_I_seq, decreasing = FALSE)
  # kparam = sort(kparam, decreasing = TRUE)

  # Combination of hyper-parameters
  # params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq[order(lambda_I_seq, decreasing = TRUE)], kparam = kparam)
  # params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq, kparam = kparam)
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            rmsvm_fit = rmlapsvm(x = x, y = y, ux = ux, gamma = gamma,
                                                 lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                                 kernel = kernel, kparam = kparam, scale = scale,
                                                 adjacency_k = adjacency_k, weightType = weightType,
                                                 normalized = normalized, ...)
                          })

                          if (!inherits(error, "try-error")){
                            pred_val = predict.rmlapsvm(rmsvm_fit, newx = valid_x)

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val$class) / length(valid_y)
                              err = 1 - acc
                            } else {
                              # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                          } else {
                            rmsvm_fit = NULL
                            err = Inf
                          }

                          return(list(error = err, fit_model = rmsvm_fit))
                        }, mc.cores = nCores)
    # valid_err = round(sapply(fold_err, "[[", "error"), 8)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  } else {
    fold_list_l = data_split(y, nfolds = nfolds)
    if (!is.null(ux)) {
      fold_list_ul = sample(rep_len(1:nfolds, length.out = nrow(ux)))
    } else {
      fold_list_ul = NULL
    }
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      fold_l = which(fold_list_l == i)
      # fold_ul = which(fold_list_ul == i)
      fold_ul = NULL
      y_fold = y[-fold_l]
      x_fold = x[-fold_l, , drop = FALSE]
      y_valid = y[fold_l]
      x_valid = x[fold_l, , drop = FALSE]
      # ux_fold = ux[-fold_ul, , drop = FALSE]
      ux_fold = ux
      fold_err = mclapply(1:nrow(params),
                          function(j) {
                            error = try({
                              msvm_fit = rmlapsvm(x = x_fold, y = y_fold, ux = ux_fold, gamma = gamma,
                                                  lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                                  kernel = kernel, kparam = kparam, scale = scale,
                                                  adjacency_k = adjacency_k, weightType = weightType,
                                                  normalized = normalized, ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.rmlapsvm(msvm_fit, newx = x_valid)
                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {
                                # err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
                              }
                            } else {
                              msvm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = msvm_fit))
                          }, mc.cores = nCores)
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
      # model_list[[i]] = lapply(fold_err, "[[", "fit_model")
    }
    # valid_err = round(colMeans(valid_err_mat), 8)
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    # opt_param = c(lambda = opt_param$lambda, lambda_I = opt_param$lambda_I)
    opt_valid_err = min(valid_err)
  }

  out$opt_param = c(lambda = opt_param$lambda, lambda_I = opt_param$lambda_I)
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$gamma = gamma
  out$scale = scale
  if (optModel) {
    opt_model = rmlapsvm(x = x, y = y, ux = ux, gamma = gamma,
                         lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
                         kernel = kernel, kparam = kparam, scale = scale,
                         adjacency_k = adjacency_k, weightType = weightType, normalized = normalized, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "rmlapsvm"
  return(out)
}
