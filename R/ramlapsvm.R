ramlapsvm_compact = function(K, L, y, gamma = 0.5, lambda, lambda_I, epsilon = 1e-6,
                             eig_tol_D = .Machine$double.eps,
                             inv_tol = 1e-25, epsilon_D = 1e-8)
{

  out = list()
  # The labeled sample size, unlabeled sample size, the number of classes and dimension of QP problem
  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)
  # if (is(y, "numeric")) {levs = as.numeric(levs)}

  # if (sum(K) == 0) {
  #   diag(K) = 1
  # }

  n = nrow(K)
  n_l = length(y_int)
  n_u = n - n_l
  qp_dim = n_l * n_class

  code_mat = code_ramsvm(y_int)
  yyi = code_mat$yyi
  W = code_mat$W
  y_index = code_mat$y_index
  Hmatj = code_mat$Hmatj
  Lmatj = code_mat$Lmatj

  J = cbind(diag(n_l), matrix(0, n_l, n_u))
  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  # LK = fixit(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K), inv_tol)
  # inv_LK = chol2inv(chol(LK))
  # K = fixit(K, eig_tol_D)

  JK = J %*% K

  LK = diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K)
  # max_LK = max(abs(LK))
  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))
  # inv_LK = solve(LK + diag(max_LK * epsilon_I, n))
  # inv_LK = solve(LK + diag(max_LK * epsilon_I, n), t(J))
  # inv_LK = inverse(LK, epsilon = inv_tol)

  # inv_LK = solve(LK / max_LK + diag(epsilon_I, n), t(J) / max_LK)
  # inv_LK = solve(LK / max_LK + diag(epsilon_I, n), tol = inv_tol / 100) / max_LK
  # inv_LK = solve(LK + diag(max(abs(diag(LK))) * epsilon_I, n), t(J), tol = inv_tol)
  # inv_LK = solve(LK, t(J), tol = inv_tol)
  inv_LK = solve(LK, tol = inv_tol) %*% t(J)
  # inv_LK = solve(LK, t(J))
  # inv_LK = solve(LK + diag(max(abs(diag(LK))) * epsilon_I, n), tol = inv_tol)
  # inv_LK = inv_LK %*% t(J)
  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))

  Q = JK %*% inv_LK
  # Q = (Q + t(Q)) / 2
  # Q = J %*% K %*% inv_LK %*% t(J)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # Q = fixit(Q, epsilon = eig_tol_D)

  # Compute Q = K x inv_LK
  D = 0
  Amat = matrix(0, n_l * n_class, n_class - 1)
  for (j in 1:(n_class - 1)) {
    D = D + Hmatj[[j]] %*% Q %*% t(Hmatj[[j]])
    Amat[, j] = -Lmatj[[j]]
  }
  D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(D))
  # D = D / max_D
  diag(D) = diag(D) + max_D * epsilon_D

  g_temp = matrix(-1, n_l, n_class)
  g_temp[y_index] = 1 - n_class
  g = as.vector(g_temp)

  dvec = -g
  # dvec = -g / max_D

  # diag(Amat[(n_class + 1):(n_class + qp_dim), ]) = 1
  # diag(Amat[(n_class + qp_dim + 1):(n_class + 2 * qp_dim), ]) = -1
  Amat = cbind(Amat, diag(-1, n_l * n_class), diag(1, n_l * n_class))

  # (3) compute Ama

  # (4) compute bvec
  # bvec = rep(0, (2 * qp_dim + n_class))

  bvec_temp = matrix(gamma - 1, nrow = n_l, ncol = n_class)
  bvec_temp[y_index] = -gamma
  if (gamma == 0 | gamma == 1) {
    bvec_temp = bvec_temp - epsilon
  }
  # bvec = c(rep(0, qp_dim + n_class), as.vector(bvec_temp))
  bvec = c(rep(0, n_class - 1), as.vector(bvec_temp), rep(0, n_l * n_class))


  # (5) find solution by solve.QP

  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, Amat, bvec, meq = n_class - 1)
  alpha = dual$solution
  alpha[alpha < 0] = 0

  alpha_mat = matrix(alpha, nrow = n_l, ncol = n_class)
  # alpha_mat[y_index][alpha_mat[y_index] > gamma] = gamma

  cmat = matrix(0, n, n_class - 1)
  for (k in 1:(n_class - 1)) {
    cmat[, k] = inv_LK %*% t(Hmatj[[k]]) %*% alpha
  }

  # find b vector using LP
  Kcmat = (JK %*% cmat) %*% W

  upper_mat = matrix(1 - gamma, nrow = nrow(alpha_mat), ncol = n_class)
  upper_mat[y_index] = gamma
  logic = (alpha_mat > 0) & (alpha_mat < upper_mat)
  if (all(colSums(logic) > 0)) {
    W_c0mat = -Kcmat - 1
    W_c0mat[y_index] = (n_class - 1) - Kcmat[y_index]
    W_c0vec = colMeans(W_c0mat)
  } else {
    alp_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
    alp_temp[y_index] = gamma

    alp = c(as.vector(alp_temp), rep(0, 2 * (n_class - 1)))


    # constraint matrix and vector
    # Alp1 = c(rep(0, qp_dim), rep(c(1, -1), n_class - 1))
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


    # print(dim(Alp))
    # print(length(blp))

    ############################################################################

    # constraint directions
    const_dir = rep(">=", qp_dim)
    # const_dir[1] = "="
    cposneg = lp("min", objective.in = alp, const.mat = Alp, const.dir = const_dir,const.rhs = blp)$solution[(qp_dim + 1):(qp_dim + 2 * (n_class - 1))]
    c0vec = rep(0, n_class - 1)
    for(j in 1:(n_class - 1)) {
      c0vec[j] = cposneg[(2 * j - 1)] - cposneg[(2 * j)]
    }
    W_c0vec = drop(t(c0vec) %*% W)
  }

  # compute the fitted values
  fit = (matrix(W_c0vec, nrow = n_l, ncol = n_class, byrow = T) + Kcmat)
  fit_class = levs[apply(fit, 1, which.max)]
  if (attr(levs, "type") == "factor") {fit_class = factor(fit_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {fit_class = as.numeric(fit_class)}
  if (attr(levs, "type") == "integer") {fit_class = as.integer(fit_class)}

  # Return the output
  out$alpha = alpha_mat
  out$cmat = cmat
  out$W_c0vec = W_c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n_l = n_l
  out$n_u = n_u
  out$n_class = n_class
  out$levels = levs
  return(out)
}


predict.ramlapsvm_compact = function(object, newK = NULL) {

  cmat = object$cmat
  W_c0vec = object$W_c0vec
  n_class = object$n_class
  levs = object$levels

  W = XI_gen(n_class)

  pred_y = matrix(W_c0vec, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  # return(list(class = pred_y, inner_prod = inner_prod))
  return(list(class = pred_class, pred_value = pred_y))
}



ramlapsvm = function(x = NULL, y, ux = NULL, gamma = 0.5, lambda, lambda_I, kernel, kparam,
                  weight = NULL, weightType = "Binary", scale = FALSE, normalized = TRUE, adjacency_k = 6, epsilon = 1e-6,
                  eig_tol_D = .Machine$double.eps, inv_tol = 1e-25, epsilon_D = 1e-8)
{
  out = list()
  n_l = NROW(x)
  n_u = NROW(ux)
  rx = rbind(x, ux)
  # X = X - matrix(colMeans(X), byrow = TRUE, nrow(X), ncol(X))
  n = n_l + n_u

  p = ncol(rx)
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

  # K = K + diag(1e-8, n)
  # K_temp = kernelMatrix(x, x, kernel = kernel, kparam = kparam) + 1

  # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # graph = W

  graph = make_knn_graph_mat(rx, k = adjacency_k)
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)

  solutions = ramlapsvm_compact(K = K, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, epsilon = epsilon,
                             eig_tol_D = eig_tol_D, inv_tol = inv_tol, epsilon_D = epsilon_D)

  out$x = x
  out$ux = ux
  out$y = y
  out$n_class = solutions$n_class
  out$levels = solutions$levels
  out$gamma = gamma
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$W_c0vec = solutions$W_c0vec
  out$epsilon = epsilon
  out$eig_tol_D = eig_tol_D
  out$inv_tol = inv_tol
  out$epsilon_D = epsilon_D
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "ramlapsvm"
  return(out)
}


predict.ramlapsvm = function(object, newx = NULL, newK = NULL, ...) {

  if (is.null(newK)) {
    newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = object$cmat
  W_c0vec = object$W_c0vec
  n_class = object$n_class
  levs = object$levels

  W = XI_gen(n_class)

  pred_y = matrix(W_c0vec, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% cmat) %*% W)
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}


cv.ramlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 10,
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
                            msvm_fit = ramlapsvm(x = x, y = y, ux = ux, gamma = gamma, lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                                 kernel = kernel, kparam = kparam, scale = scale,
                                                 adjacency_k = adjacency_k, weightType = weightType,
                                                 normalized = normalized, ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.ramlapsvm(msvm_fit, newx = valid_x)

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val$class) / length(valid_y)
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
    # valid_err = round(sapply(fold_err, "[[", "error"), 8)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    # opt_param = c(lambda = opt_param$lambda, lambda_I = opt_param$lambda_I)
    opt_valid_err = min(valid_err)
  } else {
  #   # set.seed(y[1])
    fold_list_l = data_split(y, nfolds = nfolds)
    # fold_list_ul = sample(1:nfolds, size = nrow(ux), prob = rep(1 / nfolds, nfolds), replace = TRUE)
    if (!is.null(ux)) {
      fold_list_ul = sample(rep_len(1:nfolds, length.out = nrow(ux)))
    } else {
      fold_list_ul = NULL
    }
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params), dimnames = list(paste0("Fold", 1:nfolds)))
  #   model_list = vector("list", nfolds)

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
  #     # fold = fold_list[[i]]
      fold_l = which(fold_list_l == i)
      # fold_ul = which(fold_list_ul == i)
      fold_ul = NULL
      y_fold = y[-fold_l]
      x_fold = x[-fold_l, , drop = FALSE]
      y_valid = y[fold_l]
      x_valid = x[fold_l, , drop = FALSE]
      # ux_fold = ux[-fold_ul, , drop = FALSE]
      ux_fold = ux

      #  Parallel computation on the combination of hyper-parameters
      fold_err = mclapply(1:nrow(params),
                          function(j) {
                            error = try({
                              msvm_fit = ramlapsvm(x = x_fold, y = y_fold, ux = ux_fold, gamma = gamma,
                                                   lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                                   kernel = kernel, kparam = kparam, scale = scale,
                                                   adjacency_k = adjacency_k, weightType = weightType,
                                                   normalized = normalized, ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.ramlapsvm(msvm_fit, newx = x_valid)
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
  # out$fold_models = lapply(model_list, "[[", opt_ind)
  # out$fold_ind = fold_list
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$gamma = gamma
  out$scale = scale
  if (optModel) {
    opt_model = ramlapsvm(x = x, y = y, ux = ux, gamma = gamma,
                          lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
                          kernel = kernel, kparam = kparam, scale = scale,
                          adjacency_k = adjacency_k, weightType = weightType, normalized = normalized, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "ramlapsvm"
  return(out)
}

# dyn.load("../src/alpha_update.dll")
# ramlapsvm_compact_old = function(K, L, y, gamma = 0.5, lambda, lambda_I, weight = NULL, epsilon = 1e-4 * length(y) * length(unique(y)), maxiter = 300)
# {
#
#   out = list()
#   n_class = length(unique(y))
#   n_l = length(y)
#   n = nrow(K)
#   n_u = n - n_l
#
#   inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
#   Q = (n_l * lambda) * (K %*% inv_LK)[1:n_l, 1:n_l] + 1
#
#   warm = matrix(data = 0.0, nrow = n_l, ncol = n_class)
#   if (is.null(weight)) {weight = numeric(n_l) + 1.0}
#
#   #------------------------------------------------------------------#
#   # Create k-vertex simplex.                                         #
#   #------------------------------------------------------------------#
#   my = t(XI_gen(k = n_class))
#
#   yyi = Y_matrix_gen(k = n_class, nobs = n_l, y = y)
#
#   alpha_ij = warm
#   alpha_yi = numeric(n_l)
#
#
#   erci = as.double(-rep(1, ncol(Q)) / n_l / lambda)
#
#   aa = .C("alpha_update",
#             as.vector(alpha_ij),
#             as.vector(alpha_yi),
#             as.vector(my),
#             as.vector(yyi),
#             as.vector(Q),
#             as.double(lambda),
#             as.vector(weight),
#             as.integer(n_l),
#             as.double(n_l),
#             as.integer(n_class),
#             as.double(n_class),
#             as.vector(erci),
#             as.double(gamma),
#             as.vector(as.integer(y)),
#             as.double(epsilon),
#             outalpha_ij = as.vector(numeric(n_l * n_class)),
#             maxiter = as.integer(maxiter), PACKAGE = "SMLapSVM")
#
#   warm = matrix(data = aa$outalpha_ij, nrow = n_l, ncol = n_class)
#
#   beta = beta_kernel(y = y,
#                      k = n_class,
#                      my = my,
#                      warm = warm,
#                      lambda = lambda,
#                      inv_LK = inv_LK)
#   # drop(crossprod(beta[[1]][, 1], X))
#
#   # tt = beta_linear(x = x, y = y_train, k = k, my = my, warm = warm, lambda = templambda)
#
#   # beta0 = matrix(find_theta2(y, K, gamma = gamma, cmat = beta$beta, lambda = lambda), nrow = 1)
#
#   betaout = beta$beta
#   beta0out = beta$beta0
#   # beta0out[[count]] = beta0
#
#   out$y = y
#   out$n_class = n_class
#   out$gamma = gamma
#   out$weight = weight
#   out$lambda = lambda
#   out$lambda_I = lambda_I
#   out$beta = betaout
#   out$beta0 = beta0out
#   out$epsilon = epsilon
#   out$warm = warm
#   class(out) = "ramlapsvm_compact"
#   return(out)
# }

# predict.ramlapsvm_compact = function(object, newK = NULL) {
#
#   beta = object$beta
#   beta0 = object$beta0
#
#   temp_pred_y = predict_kernel(K_test = newK,
#                                beta = beta,
#                                beta0 = beta0,
#                                k = object$n_class)
#   inner_prod = temp_pred_y$inner_prod
#   temp_pred_y = temp_pred_y$class
#   pred_y = temp_pred_y
#
#   return(list(class = pred_y, inner_prod = inner_prod))
# }




