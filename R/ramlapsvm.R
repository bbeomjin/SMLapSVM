# dyn.load("../src/alpha_update.dll")
ramlapsvm = function(x = NULL, y, ux = NULL, gamma = 0.5, lambda, lambda_I, kernel, kparam,
                  weight = NULL, weightType = "Binary", scale = FALSE, normalized = TRUE, adjacency_k = 6,
                  epsilon = 1e-4 * length(y) * length(unique(y)), warm = NULL,
                  maxiter = 300) {

  n_l = NROW(x)
  n_u = NROW(ux)
  rx = rbind(x, ux)
  # X = X - matrix(colMeans(X), byrow = TRUE, nrow(X), ncol(X))
  n = n_l + n_u
  n_class = length(unique(y))

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

  K = NULL
  if (is.null(K)) {
    K = kernelMat(rx, rx, kernel = kernel, kparam = kparam)
    # K = kernelMatrix(rbfdot(sigma = kparam), x, x) + 1
  } else {
    K = K
  }
  # K = K + diag(1e-8, n)
  # K_temp = kernelMat(x, x, kernel = kernel, kparam = kparam) + 1

  # graph = make_knn_graph_mat(rx, k = adjacency_k)
  W = RSSL:::adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  graph = W
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType)

  # W = RSSL:::adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # d = rowSums(W)
  # L = diag(d) - W
  # if (normalized) {
  #   L = diag(1 / sqrt(d)) %*% L %*% diag(1 / sqrt(d))
  # }

  #------------------------------------------------------------------#
  # Make Q matrix.                                                   #
  #------------------------------------------------------------------#
  inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))

  # as.vector(((temp_L %*% solve(K + diag(1e-10, n))) %*% K[, 1]))
  # temp_L[1, ]
  # temp_L[, 1]

  # a = matrix(0, nrow = 100, ncol = 100)
  # a[upper.tri(a)] = rnorm(n = length(a[upper.tri(a)]))
  # a = t(a) + a
  # b = matrix(0, nrow = 100, ncol = 100)
  # b[upper.tri(b)] = rnorm(n = length(b[upper.tri(b)]))
  # b = t(b) + b
  # a %*% b == t(b %*% a)

  Q = (n_l * lambda) * (K %*% inv_LK)[1:n_l, 1:n_l] + 1

  #------------------------------------------------------------------#
  # Convert labels to integers.                                      #
  #------------------------------------------------------------------#

  #------------------------------------------------------------------#
  # Calculate kernel                                                 #
  #------------------------------------------------------------------#

  warm = matrix(data = 0.0, nrow = n_l, ncol = n_class)
  if (is.null(weight)) {weight = numeric(n_l) + 1.0}

  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  my = t(XI_gen(k = n_class))

  yyi = Y_matrix_gen(k = n_class,
                     nobs = n_l,
                     y = y)

  alpha_ij = warm
  alpha_yi = numeric(n_l)

  # erci = -diag(Q) / 2 / nobsdouble / templambda
  erci = -as.double(rep(1, ncol(Q)) / n_l / lambda)

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

  out = list()
  out$x = x
  out$ux = ux
  out$y = y
  out$n_class = n_class
  out$gamma = gamma
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$kparam = kparam
  out$beta = betaout
  out$beta0 = beta0out
  out$epsilon = epsilon
  out$warm = warm
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  class(out) = "ramlapsvm"
  return(out)
}


predict.ramlapsvm = function(object, newx = NULL, newK = NULL, ...) {

  if (is.null(newK)) {
    newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

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



Kfold_ramlapsvm = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 10,
                         lambda_seq = 2^{seq(-10, 15, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                         scale = FALSE, adjacency_k = 6, weightType = "Heatmap", normalized = TRUE,
                         gamma = 0.5, kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                         criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
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
                          msvm_fit = ramlapsvm(x = x, y = y, ux = ux, gamma = gamma, lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                             kernel = kernel, kparam = params$kparam[j], scale = scale, adjacency_k = adjacency_k,
                                             weightType = weightType, normalized = normalized, ...)
                          pred_val = predict.ramlapsvm(msvm_fit, newx = valid_x)

                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val[[1]][[1]]) / length(valid_y)
                            err = 1 - acc
                          } else {
                            err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                          }
                          return(list(error = err, fit_model = msvm_fit))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }
  # } else {
  #   # set.seed(y[1])
  #   # fold_list = createFolds(y, k = nfolds, list = TRUE)
  #   fold_list = data_split(y, nfolds, k = k)
  #   valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
  #   model_list = vector("list", nfolds)
  #
  #   for (i in 1:nfolds) {
  #     cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
  #     # fold = fold_list[[i]]
  #     fold = which(fold_list == i)
  #     y_fold = y[-fold]
  #     x_fold = x[-fold, , drop = FALSE]
  #     y_valid = y[fold]
  #     x_valid = x[fold, , drop = FALSE]
  #
  #     #  Parallel computation on the combination of hyper-parameters
  #     fold_err = mclapply(1:nrow(params),
  #                         function(j) {
  #                           msvm_fit = SRAMSVM_solve(x = x_fold, y = y_fold, gamma = gamma,
  #                                                    lambda = params$lambda[j], kernel = kernel,
  #                                                    kparam = params$kparam[j], ...)
  #                           pred_val = predict(msvm_fit, newx = x_valid)
  #
  #                           if (criterion == "0-1") {
  #                             acc = sum(y_valid == pred_val[[1]][[1]]) / length(y_valid)
  #                             err = 1 - acc
  #                           } else {
  #                             err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
  #                           }
  #                           return(list(error = err, fit_model = msvm_fit))
  #                         }, mc.cores = nCores)
  #     valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
  #     model_list[[i]] = lapply(fold_err, "[[", "fit_model")
  #   }
  #   valid_err = colMeans(valid_err_mat, na.rm = TRUE)
  #   opt_ind = max(which(valid_err == min(valid_err)))
  #   opt_param = params[opt_ind, ]
  #   opt_valid_err = min(valid_err)
  # }

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
  out$gamma = gamma
  out$scale = scale
  if (optModel) {
    opt_model = ramlapsvm(x = x, y = y, ux = ux, gamma = gamma, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
                        kernel = kernel, kparam = opt_param$kparam, scale = scale, adjacency_k = adjacency_k,
                        weightType = weightType, normalized = normalized, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "ramlapsvm"
  return(out)
}


# For parallel computing in windows  --------------------------------------

# mclapply.hack = function(..., nCores) {
#   ## Create a cluster
#   size.of.list = length(list(...)[[1]])
#   cl <- makeCluster(nCores)
#   ## Find out the names of the loaded packages
#   loaded.package.names = c(
#     ## Base packages
#     sessionInfo()$basePkgs,
#     ## Additional packages
#     names( sessionInfo()$otherPkgs ))
#   tryCatch({
#     ## Copy over all of the objects within scope to
#     ## all clusters.
#     this.env = environment()
#     while (identical(this.env, globalenv()) == FALSE) {
#       clusterExport(cl,
#                     ls(all.names = TRUE, env = this.env),
#                     envir = this.env)
#       this.env = parent.env(environment())
#     }
#     clusterExport(cl,
#                   ls(all.names = TRUE, env = globalenv()),
#                   envir = globalenv())
#
#     ## Load the libraries on all the clusters
#     ## N.B. length(cl) returns the number of clusters
#     parLapply(cl, 1:length(cl), function(xx){
#       lapply(loaded.package.names, function(yy) {
#         require(yy, character.only = TRUE)})
#     })
#
#     ## Run the lapply in parallel
#     return(parLapply(cl, ...) )
#   }, finally = {
#     ## Stop the cluster
#     stopCluster(cl)
#   })
# }
#
#
# Kfold_mlapsvm_angle_win = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 10,
#                          lambda_seq = 2^{seq(-10, 15, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
#                          gamma = 0.5, scale = FALSE, adjacency_k = 6, weightType = "Heatmap", normalized = TRUE,
#                          kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
#                           criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
# {
#   call = match.call()
#   kernel = match.arg(kernel)
#   criterion = match.arg(criterion)
#
#   # if (scale) {
#   #   x = scale(x)
#   #   if (!is.null(valid_x)) {
#   #     means = attr(x, "scaled:center")
#   #     stds = attr(x, "scaled:scale")
#   #     valid_x = (valid_x - matrix(means, NROW(x), NCOL(x), byrow = TRUE)) / matrix(stds, NROW(x), NCOL(x), byrow = TRUE)
#   #   }
#   # }
#
#   if (!is.numeric(lambda_seq)) {
#     lambda_seq = as.numeric(lambda_seq)
#   }
#
#   if (!is.numeric(lambda_I_seq)) {
#     lambda_I_seq = as.numeric(lambda_I_seq)
#   }
#
#   if (!is.numeric(kparam)) {
#     kparam = as.numeric(kparam)
#   }
#
#   # The number of classes
#   k = length(unique(y))
#
#   # Combination of hyper-parameters
#   params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq, kparam = kparam)
#
#   if (!is.null(valid_x) & !is.null(valid_y)) {
#     model_list = vector("list", 1)
#     fold_list = NULL
#
#     #  Parallel computation on the combination of hyper-parameters
#     fold_err = mclapply.hack(1:nrow(params),
#                         function(j) {
#                           msvm_fit = mlapsvm_angle(x = x, y = y, ux = ux, gamma = gamma, lambda = params$lambda[j], lambda_I = params$lambda_I[j],
#                                              kernel = kernel, kparam = params$kparam[j], scale = scale, adjacency_k = adjacency_k,
#                                              weightType = weightType, normalized = normalized, ...)
#                           pred_val = predict.mlapsvm_angle(msvm_fit, newx = valid_x)
#
#                           if (criterion == "0-1") {
#                             acc = sum(valid_y == pred_val[[1]][[1]]) / length(valid_y)
#                             err = 1 - acc
#                           } else {
#                             err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
#                           }
#                           return(list(error = err, fit_model = msvm_fit))
#                         }, nCores = nCores)
#     valid_err = sapply(fold_err, "[[", "error")
#     model_list[[1]] = lapply(fold_err, "[[", "fit_model")
#     opt_ind = max(which(valid_err == min(valid_err)))
#     opt_param = params[opt_ind, ]
#     opt_valid_err = min(valid_err)
#   }
#   # } else {
#   #   # set.seed(y[1])
#   #   # fold_list = createFolds(y, k = nfolds, list = TRUE)
#   #   fold_list = data_split(y, nfolds, k = k)
#   #   valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params))
#   #   model_list = vector("list", nfolds)
#   #
#   #   for (i in 1:nfolds) {
#   #     cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
#   #     # fold = fold_list[[i]]
#   #     fold = which(fold_list == i)
#   #     y_fold = y[-fold]
#   #     x_fold = x[-fold, , drop = FALSE]
#   #     y_valid = y[fold]
#   #     x_valid = x[fold, , drop = FALSE]
#   #
#   #     #  Parallel computation on the combination of hyper-parameters
#   #     fold_err = mclapply(1:nrow(params),
#   #                         function(j) {
#   #                           msvm_fit = SRAMSVM_solve(x = x_fold, y = y_fold, gamma = gamma,
#   #                                                    lambda = params$lambda[j], kernel = kernel,
#   #                                                    kparam = params$kparam[j], ...)
#   #                           pred_val = predict(msvm_fit, newx = x_valid)
#   #
#   #                           if (criterion == "0-1") {
#   #                             acc = sum(y_valid == pred_val[[1]][[1]]) / length(y_valid)
#   #                             err = 1 - acc
#   #                           } else {
#   #                             err = ramsvm_hinge(y_valid, pred_val$inner_prod, k = k, gamma = gamma)
#   #                           }
#   #                           return(list(error = err, fit_model = msvm_fit))
#   #                         }, mc.cores = nCores)
#   #     valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
#   #     model_list[[i]] = lapply(fold_err, "[[", "fit_model")
#   #   }
#   #   valid_err = colMeans(valid_err_mat, na.rm = TRUE)
#   #   opt_ind = max(which(valid_err == min(valid_err)))
#   #   opt_param = params[opt_ind, ]
#   #   opt_valid_err = min(valid_err)
#   # }
#
#   out = list()
#   out$opt_param = opt_param
#   out$opt_valid_err = opt_valid_err
#   out$opt_ind = opt_ind
#   out$valid_err = valid_err
#   out$fold_models = lapply(model_list, "[[", opt_ind)
#   out$fold_ind = fold_list
#   out$x = x
#   out$y = y
#   out$valid_x = valid_x
#   out$valid_y = valid_y
#   out$kernel = kernel
#   out$gamma = gamma
#   out$scale = scale
#   if (optModel) {
#     opt_model = mlapsvm_angle(x = x, y = y, ux = ux, gamma = gamma, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
#                         kernel = kernel, kparam = opt_param$kparam, scale = scale, adjacency_k = adjacency_k,
#                         weightType = weightType, normalized = normalized, ...)
#     out$opt_model = opt_model
#   }
#   out$call = call
#   class(out) = "mlapsvm"
#   return(out)
# }
