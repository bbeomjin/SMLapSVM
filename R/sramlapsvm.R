sramlapsvm = function(x = NULL, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                      lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                      lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                      gamma = 0.5, adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                      kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                      scale = FALSE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, verbose = 0, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep.sramlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                               lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq,
                               theta = NULL, fold_theta = NULL,
                               gamma = gamma, adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                               kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  cat("Fit theta-step \n")
  thetastep_fit = thetastep.sramlapsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.sramlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                                   lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq,
                                   theta = thetastep_fit$opt_theta, fold_theta = thetastep_fit$opt_fold_theta,
                                   gamma = gamma, adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                                   kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)

  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
    cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
    cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
  }

  out$opt_theta = thetastep_fit$opt_theta
  out$cstep_inform = cstep_fit
  out$thetastep_inform = thetastep_fit
  out$opt_param = opt_cstep_fit$opt_param
  out$opt_model = opt_cstep_fit$opt_model
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$valid_err = opt_cstep_fit$valid_err
  out$type = type
  out$call = call
  class(out) = "sramlapsvm"
  return(out)
}

predict.sramlapsvm = function(object, newx = NULL, newK = NULL)
{
  model = object$opt_model
  beta = model$beta
  beta0 = model$beta0
  levs = model$levels

  # if (object$scale) {
  #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  # }

  if (is.null(newK)) {
    new_anova_K = make_anovaKernel(newx, rbind(object$cstep_inform$x, object$cstep_inform$ux),
                                   kernel = object$cstep_inform$kernel, kparam = object$cstep_inform$kparam)
    newK = combine_kernel(new_anova_K, theta = object$opt_theta)
    # newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  W = XI_gen(model$n_class)

  W_beta0 = drop(t(beta0) %*% W)

  pred_y = matrix(W_beta0, nrow = nrow(newK), ncol = model$n_class, byrow = T) + ((newK %*% beta) %*% W)
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}


cstep.sramlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                            lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                            theta = NULL, fold_theta = NULL,
                            kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                            scale = FALSE, adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                            criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  p = ncol(x)

  lambda_seq = as.numeric(lambda_seq)
  lambda_I_seq = as.numeric(lambda_I_seq)
  kparam = as.numeric(kparam)

  # 추후, 커널에 맞게 theta의 길이 조절
  if (is.null(theta)) {
    theta = rep(1, p)
  }

  if (is.null(fold_theta)) {
    fold_theta = rep(list(rep(1, p)), nfolds)
  }

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  lambda_I_seq = sort(lambda_I_seq, decreasing = FALSE)
  # kparam = sort(kparam, decreasing = FALSE)

  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq)

  n_l = NROW(x)
  n_u = NROW(ux)
  n = n_l + n_u
  rx = rbind(x, ux)

  center = rep(0, p)
  scaled = rep(1, p)

  if (scale) {
    rx = scale(rx)
    center = attr(rx, "scaled:center")
    scaled = attr(rx, "scaled:scale")
    x = (x - matrix(center, nrow = n_l, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_l, ncol = p, byrow = TRUE)
    ux = (ux - matrix(center, nrow = n_u, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_u, ncol = p, byrow = TRUE)
  }

  graph = make_knn_graph_mat(rx, k = adjacency_k)
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = kparam)
    # K = combine_kernel(anova_kernel = anova_K, theta = theta)

    # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
    # graph = W

    # L = fixit(L, epsilon = 0)

    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = kparam)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)
    #  Parallel computation on the combination of hyper-parameters

    # fit = angle_lapsvm_core(K = K, L = L, y = y, lambda = params$lambda[j], lambda_I = params$lambda_I[j], maxiter = 1e+4)

    # predict.angle_lapsvm_core(fit, newK = K[1:30, ])[[2]]
    # table(predict.angle_lapsvm_core(fit, newK = K[1:30, ])[[1]], y)

    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            msvm_fit = sramlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                          lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                            # msvm_fit = angle_lapsvm_core(K = K, L = L, y = y, lambda = params$lambda[j], lambda_I = params$lambda_I[j], gamma = gamma)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.ramlapsvm_compact(msvm_fit, newK = valid_K)

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
    valid_err = round(sapply(fold_err, "[[", "error"), 8)
    # valid_err = sapply(fold_err, "[[", "error")
    # model_list[[1]] = lapply(fold_err, "[[", "fit_model")
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

    fold_list = list(labeled = fold_list_l, unlabeled = fold_list_ul)

    valid_err_mat = matrix(NA, nrow = nfolds, ncol = nrow(params), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      #     # fold = fold_list[[i]]
      fold_l = which(fold_list_l == i)
      # fold_ul = which(fold_list_ul == i)
      fold_ul = NULL
      y_train = y[-fold_l]
      x_train = x[-fold_l, , drop = FALSE]
      y_valid = y[fold_l]
      x_valid = x[fold_l, , drop = FALSE]
      # ux_train = ux[-fold_ul, , drop = FALSE]
      ux_train = ux
      rx_train = rbind(x_train, ux_train)

      theta_train = fold_theta[[i]]

      subanova_K = make_anovaKernel(rx_train, rx_train, kernel, kparam)
      # subK = combine_kernel(subanova_K, theta)

      subanova_K_valid = make_anovaKernel(x_valid, rx_train, kernel, kparam)
      subK_valid = combine_kernel(subanova_K_valid, theta_train)


      graph_train = make_knn_graph_mat(rx_train, k = adjacency_k)
      L_train = make_L_mat(rx_train, kernel = kernel, kparam = kparam, graph = graph_train, weightType = weightType, normalized = normalized)

      fold_err = mclapply(1:nrow(params),
                          function(j) {
                            error = try({
                              msvm_fit = sramlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta_train, y = y_train, gamma = gamma,
                                                            lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.ramlapsvm_compact(msvm_fit, newK = subK_valid)

                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {

                              }
                            } else {
                              msvm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = msvm_fit))
                          }, mc.cores = nCores)
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
    }
    valid_err = round(colMeans(valid_err_mat), 8)
    # valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }
  out$opt_param = c(lambda = opt_param$lambda, lambda_I = opt_param$lambda_I)
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$ux = ux
  out$y = y
  out$L = L
  out$adjacency_k = adjacency_k
  out$normalized = normalized
  out$weightType = weightType
  out$theta = theta
  out$fold_theta = fold_theta
  out$gamma = gamma
  out$valid_x = valid_x
  out$valid_y = valid_y
  # out$anova_K = anova_K
  # out$K = K
  # out$valid_anova_K = valid_anova_K
  # out$valid_K = valid_K
  out$kernel = kernel
  out$kparam = kparam
  out$scale = scale
  out$nfolds = nfolds
  out$fold_list = fold_list
  out$criterion = criterion
  if (optModel) {
    anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = kparam)
    opt_model = sramlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                   lambda = out$opt_param["lambda"], lambda_I = out$opt_param["lambda_I"], ...)
    # opt_model = angle_lapsvm_core(K = K, L = L, y = y, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, gamma = gamma)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "sramlapsvm"
  return(out)
}

thetastep.sramlapsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                                isCombined = TRUE, optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  lambda_I = object$opt_param["lambda_I"]
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$kparam
  gamma = object$gamma
  x = object$x
  y = object$y
  init_theta = object$theta
  init_fold_theta = object$fold_theta
  ux = object$ux
  rx = rbind(x, ux)

  valid_x = object$valid_x
  valid_y = object$valid_y
  adjacency_k = object$adjacency_k
  normalized = object$normalized
  weightType = object$weightType
  L = object$L
  nfolds = object$nfolds
  fold_list = object$fold_list
  # anova_K = object$anova_K
  # K = object$K
  # valid_anova_K = object$valid_anova_K

  anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = kparam)

  if (is.null(object$opt_model)) {
    opt_model = sramlapsvm_compact(anova_K = anova_K, L = L, theta = init_theta, y = y, lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
  } else {
    opt_model = object$opt_model
  }

  if (!is.null(valid_x) & !is.null(valid_y)) {

    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = kparam)

    init_model = opt_model

    fold_err = mclapply(1:length(lambda_theta_seq),
                        function(j) {
                          error = try({
                            theta = find_theta.sramlapsvm(y = y, gamma = gamma, anova_kernel = anova_K, L = L,
                                                          cmat = init_model$beta, c0vec = init_model$beta0,
                                                          lambda = lambda, lambda_I = lambda_I,
                                                          lambda_theta = lambda_theta_seq[j], ...)

                            if (isCombined) {
                              # subK = combine_kernel(anova_K, theta)
                              init_model = sramlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                              lambda = lambda, lambda_I = lambda_I, ...)
                              # init_model = angle_lapsvm_core(K = subK, L = L, y = y, lambda = lambda, lambda_I = lambda_I, gamma = gamma)
                            }
                          })

                          if (!inherits(error, "try-error")) {
                            valid_subK = combine_kernel(valid_anova_K, theta)
                            pred_val = predict.ramlapsvm_compact(init_model, newK = valid_subK)

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val$class) / length(valid_y)
                              err = 1 - acc
                            } else {
                              # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                          } else {
                            err = Inf
                            theta = rep(0, anova_K$numK)
                          }

                          return(list(error = err, theta = theta))
                        }, mc.cores = nCores)
    valid_err = round(sapply(fold_err, "[[", "error"), 8)
    # valid_err = sapply(fold_err, "[[", "error")
    theta_seq = sapply(fold_err, "[[", "theta")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_theta = theta_seq[, opt_ind]
    opt_valid_err = min(valid_err)
    opt_fold_theta = init_fold_theta
  } else {
    fold_list_l = fold_list$labeled
    fold_list_ul = fold_list$unlabeled

    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq), dimnames = list(paste0("Fold", 1:nfolds)))
    fold_theta = vector("list", nfolds)
    names(fold_theta) = paste0("Fold", 1:nfolds)
    # theta_seq_list = mclapply(1:length(lambda_theta_seq),
    #                           function(j) {
    #                             error = try({
    #                               theta = find_theta.sramlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = opt_model$beta, c0vec = opt_model$beta0,
    #                                                             gamma = gamma, lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
    #                             })
    #                             if (inherits(error, "try-error")) {
    #                               theta = rep(0, anova_K$numK)
    #                             }
    #                             return(theta)
    #                           }, mc.cores = nCores)

    for (i in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      #     # fold = fold_list[[i]]
      fold_l = which(fold_list_l == i)
      # fold_ul = which(fold_list_ul == i)
      fold_ul = NULL
      y_train = y[-fold_l]
      x_train = x[-fold_l, , drop = FALSE]
      y_valid = y[fold_l]
      x_valid = x[fold_l, , drop = FALSE]
      # ux_train = ux[-fold_ul, , drop = FALSE]
      ux_train = ux
      rx_train = rbind(x_train, ux_train)

      init_theta_train = init_fold_theta[[i]]

      subanova_K = make_anovaKernel(rx_train, rx_train, kernel, kparam)
      # subK = combine_kernel(subanova_K, theta)

      subanova_K_valid = make_anovaKernel(x_valid, rx_train, kernel, kparam)
      # subK_valid = combine_kernel(subanova_K_valid, theta)

      graph_train = make_knn_graph_mat(rx_train, k = adjacency_k)
      L_train = make_L_mat(rx_train, kernel = kernel, kparam = kparam, graph = graph_train, weightType = weightType, normalized = normalized)

      init_model = sramlapsvm_compact(anova_K = subanova_K, L = L_train, theta = init_theta_train, y = y_train,
                                      lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
      cmat = init_model$beta
      c0vec = init_model$beta0

      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = find_theta.sramlapsvm(y = y_train, anova_kernel = subanova_K, L = L_train, cmat = cmat, c0vec = c0vec,
                                                            gamma = gamma, lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)

                              if (isCombined) {
                                init_model = sramlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta, y = y_train,
                                                                lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
                              }
                            })

                            # error = try({
                            #   theta = theta_seq_list[[j]]
                            #   if (isCombined) {
                            #     init_model = sramlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta, y = y_train,
                            #                                     lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
                            #   }
                            # })

                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.ramlapsvm_compact(init_model, newK = subK_valid)

                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {

                              }
                            } else {
                              err = Inf
                              theta = rep(0, subanova_K$numK)
                            }
                            return(list(theta = theta, error = err))
                          }, mc.cores = nCores)
      fold_theta[[i]] = do.call(cbind, lapply(fold_err, "[[", "theta"))
      valid_err_mat[i, ] = sapply(fold_err, "[[", "error")
    }
    valid_err = round(colMeans(valid_err_mat), 8)
    # valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)

    opt_fold_theta = lapply(fold_theta, function(x) return(x[, opt_ind]))

    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = find_theta.sramlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = opt_model$beta, c0vec = opt_model$beta0,
                                                                gamma = gamma, lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
                                })
                                if (inherits(error, "try-error")) {
                                  theta = rep(0, anova_K$numK)
                                }
                                return(theta)
                              }, mc.cores = nCores)

    theta_seq = do.call(cbind, theta_seq_list)
    opt_theta = theta_seq[, opt_ind]
  }

  out$opt_lambda_theta = opt_lambda_theta
  out$opt_ind = opt_ind
  out$opt_theta = opt_theta
  out$theta_seq = theta_seq
  out$opt_valid_err = opt_valid_err
  out$valid_err = valid_err
  out$opt_fold_theta = opt_fold_theta
  if (optModel) {
    # subK = combine_kernel(anova_K, opt_theta)
    opt_model = sramlapsvm_compact(anova_K = anova_K, L = L, theta = opt_theta, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
    out$opt_model = opt_model
  }
  class(out) = "sramlapsvm"
  return(out)
}


find_theta.sramlapsvm = function(y, anova_kernel, L, cmat, c0vec, gamma, lambda, lambda_I, lambda_theta = 1,
                                 eig_tol_D = .Machine$double.eps,
                                 eig_tol_I = .Machine$double.eps,
                                 epsilon_D = 1e-8,
                                 epsilon_I = .Machine$double.eps,
                                 inv_tol = 1e-25)
{
  if (lambda_theta <= 0) {
    theta = rep(1, anova_kernel$numK)
    return(theta)
  }

  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }

  n = NROW(cmat)

  # anova_kernel_orig = anova_kernel
  # anova_kernel$K = lapply(anova_kernel$K, function(x) {
  #   diag(x) = diag(x) + max(abs(diag(x))) * epsilon_I
  #   return(x)
  # })

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)
  n_class = length(levs)

  # Standard QP form :
  # min
  # n = NROW(cmat)
  n_l = length(y_int)
  n_u = n - n_l

  trans_Y = Y_matrix_gen(n_class, nobs = n_l, y = y_int)
  Y_code = XI_gen(n_class)
  # Y_code = Y_matrix_gen(n_class, nobs = n_class, y = 1:n_class)
  y_index = cbind(1:n_l, y_int)

  cmat = as.matrix(cmat)
  c0vec = as.vector(c0vec)

  Dmat = numeric(anova_kernel$numK)
  dvec = numeric(anova_kernel$numK)
  for (j in 1:anova_kernel$numK) {
    temp_D = 0
    temp_d = 0
    # temp_A = NULL
    for (q in 1:(n_class - 1)) {
      cvec = cmat[, q]
      # KLK_temp = anova_kernel_orig$K[[j]] %*% L %*% anova_kernel_orig$K[[j]]
      # diag(KLK_temp) = diag(KLK_temp) + max(abs(KLK_temp)) * epsilon_I

      temp_D = temp_D + n_l * lambda_I / (1 * n^2) * t(cvec) %*% anova_kernel$K[[j]] %*% L %*% anova_kernel$K[[j]] %*% cvec
      temp_d = temp_d + n_l * lambda / 2 * t(cvec) %*% anova_kernel$K[[j]] %*% cvec + n_l * lambda_theta
    }
    Dmat[j] = temp_D
    dvec[j] = temp_d
  }

  A_mat = NULL
  for (j in 1:anova_kernel$numK) {
    temp_A = NULL
    for (k in 1:n_class) {
      inner_prod_A = matrix(rowSums((anova_kernel$K[[j]][1:n_l, ] %*% cmat) * matrix(Y_code[, k], nrow = n_l, ncol = n_class - 1, byrow = TRUE)))
      temp_A = rbind(temp_A, inner_prod_A)
    }
    A_mat = cbind(A_mat, temp_A)
  }

  max_D = max(abs(Dmat))
  Dmat = c(Dmat, c(rep(0, n_l * n_class)))
  Dmat = diag(Dmat)
  # Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)
  diag(Dmat) = diag(Dmat) + max_D * epsilon_D
  # diag(Dmat) = diag(Dmat) + epsilon_D

  # Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)
  # diag(Dmat) = diag(Dmat) + 1e-8
  # Dmat = Dmat / max_D

  dvec_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  dvec_temp[y_index] = gamma
  # dvec_temp = as.vector(Y)
  # dvec_temp[dvec_temp == 1] = 0
  # dvec_temp[dvec_temp < 0] = 1
  dvec = c(dvec, as.vector(dvec_temp))
  # dvec = dvec / max_D
  # solve QP
  # diag(Dmat) = diag(Dmat) + epsilon

  m_index = matrix(1:(n_l * n_class), ncol = n_class)[y_index]
  A_mat[m_index, ] = -A_mat[m_index, ]
  A_mat = cbind(-A_mat, diag(1, n_l * n_class))
  A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)
  #    print(ncol(A.mat))

  bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * c0vec[x]))
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(rowSums(t(Y_code) * matrix(c0vec, nrow = ncol(Y_code), ncol = nrow(Y_code), byrow = TRUE)), nrow = n_l, ncol = n_class, byrow = TRUE)
  bb_j[y_index] = bb_yi
  bvec = c(as.vector(bb_j), rep(0, anova_kernel$numK + n_l * n_class), rep(-1, anova_kernel$numK))

  #    print(A.mat)
  #    print(bvec)
  theta_sol = solve.QP(Dmat, -dvec, t(A_mat), bvec, meq = 0, factorized = FALSE)$solution
  theta = theta_sol[1:anova_kernel$numK]
  # theta[theta < 1e-8] = 0
  # theta = round(theta, 6)
  # theta_sol[theta_sol < 1e-6] = 0
  #    print(beta)
  return(theta)
}


sramlapsvm_compact = function(anova_K, L, theta, y, gamma = 0.5, lambda, lambda_I, epsilon = 1e-6,
                              eig_tol_D = .Machine$double.eps,
                              eig_tol_I = .Machine$double.eps,
                              epsilon_D = 1e-8,
                              epsilon_I = .Machine$double.eps,
                              inv_tol = 1e-25)
{

  out = list()
  # The labeled sample size, unlabeled sample size, the number of classes and dimension of QP problem
  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)

  # n = nrow(anova_K$K[[1]])

  # anova_K_orig = anova_K
  # anova_K$K = lapply(anova_K$K, function(x) {
  #   diag(x) = diag(x) + max(abs(diag(x))) * epsilon_I
  #   return(x)
  # })

  K = combine_kernel(anova_K, theta = theta)
  n = nrow(K)

  if (sum(K) == 0) {
    diag(K) = 1
  }

  # n = nrow(K)
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

  # KLK = 0
  # for (i in 1:anova_K$numK) {
  #   KLK_temp = anova_K_orig$K[[i]] %*% L %*% anova_K_orig$K[[i]]
  #   diag(KLK_temp) = diag(KLK_temp) + max(abs(KLK_temp)) * epsilon_I
  #   KLK = KLK + theta[i]^2 * KLK_temp
  # }
  KLK = 0
  for (i in 1:anova_K$numK) {
    KLK = KLK + theta[i]^2 * anova_K$K[[i]] %*% L %*% anova_K$K[[i]]
  }
  # KLK = (KLK + t(KLK)) / 2

  # K = fixit(K, epsilon = eig_tol_I)
  # KLK = fixit(KLK, epsilon = eig_tol_I)

  # lambda_K = n_l * lambda * K
  # lambda_K = fixit(lambda_K, epsilon = eig_tol_I)
  # diag(lambda_K) = diag(lambda_K) + epsilon_I
  # diag(lambda_K) = diag(lambda_K) + max(abs(diag(lambda_K))) * epsilon_I

  # lambda_KLK = n_l * lambda_I / n^2 * KLK
  # lambda_KLK = fixit(lambda_KLK, epsilon = eig_tol_I)
  # lambda_KLK = fixit(lambda_KLK)
  # diag(lambda_KLK) = diag(lambda_KLK) + epsilon_I
  # diag(lambda_KLK) = diag(lambda_KLK) + max(abs(diag(lambda_KLK))) * epsilon_I

  # K_KLK = lambda_K + lambda_KLK
  K_KLK = n_l * lambda * K + n_l * lambda_I / n^2 * KLK
  # K_KLK = (K_KLK + t(K_KLK)) / 2
  K_KLK = fixit(K_KLK, epsilon = eig_tol_I)
  diag(K_KLK) = diag(K_KLK) + max(abs(diag(K_KLK))) * epsilon_I
  # diag(K_KLK) = diag(K_KLK) + epsilon_I

  JK = J %*% K

  # inv_K_KLK = solve(K_KLK, tol = inv_tol)
  # inv_K_KLK = chol2inv(chol(K_KLK))
  # inv_K_KLK = inverse(K_KLK)
  # inv_K_KLK = (inv_K_KLK + t(inv_K_KLK)) / 2
  # inv_K_KLK = tcrossprod(inv_K_KLK, JK)
  # inv_K_KLK = inv_K_KLK %*% t(JK)
  inv_K_KLK = solve(K_KLK, t(JK), tol = inv_tol)
  # inv_K_KLK = solve(K_KLK, tol = inv_tol)
  # inv_K_KLK = qr.solve(K_KLK, K %*% t(J), tol = eig_tol_I)

  # Q = JK %*% inv_K_KLK %*% t(JK)
  Q = JK %*% inv_K_KLK
  Q = (Q + t(Q)) / 2
  # Q = J %*% K %*% inv_KLK

  # Q = fixit(Q, epsilon = eig_tol_D)
  # diag(Q) = diag(Q) + epsilon_D

  # Compute Q = K x inv_LK
  D = 0
  Amat = matrix(0, n_l * n_class, n_class - 1)
  for (j in 1:(n_class - 1)) {
    D = D + Hmatj[[j]] %*% Q %*% t(Hmatj[[j]])
    Amat[, j] = -Lmatj[[j]]
  }
  # D = (D + t(D)) / 2
  # D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(diag(D)))
  # D = D / max_D
  diag(D) = diag(D) + max_D * epsilon_D
  # diag(D) = diag(D) + epsilon_D
  #################################### for test #######################################
  # alpha_mat = matrix(rnorm(n_l * n_class), n_l, n_class)
  # temp_vec = 0
  # for (i in 1:n_l) {
  #   y_temp = y[i]
  #   temp_vec = temp_vec + (W[1, y_temp] * alpha_mat[i, y_temp] * K[, i] - rowSums(matrix(W[1, -y_temp] * alpha_mat[i, -y_temp], nrow = n, ncol = n_class - 1, byrow = TRUE) * matrix(K[, i], nrow = n, ncol = n_class - 1)))
  # }
  # temp = as.vector(alpha_mat) %*% Hmatj[[1]] %*% J %*% K
  #################################### for test #######################################


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

  # for (j in 1:n_class) {
  #   for (i in 1:n_l) {
  #     flag = 0
  #     if (y[i] == j) {
  #       flag = 1
  #     }
  #     bvec[n_class + qp_dim + (j - 1) * n_l + i] = -(gamma * flag + (1 - gamma) * (1 - flag))
  #     # correction to avoid redundant constraints when gamma = 0 or 1
  #     if ((flag == 1 & gamma == 0) | (flag == 0 & gamma == 1)) {
  #       bvec[n_class + qp_dim + (j - 1) * n_l + i] = bvec[n_class + qp_dim + (j - 1) * n_l + i] - epsilon
  #     }
  #   }
  # }

  # remove one redudant constraint
  # Amat1 = Amat[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class)), ]
  # bvec1 = bvec[c(1:(n_class - 1), (n_class + 1):(2 * qp_dim + n_class))]

  # (5) find solution by solve.QP

  nonzero = find_nonzero(Amat)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(D, dvec, Amat, Aind, bvec, meq = (n_class - 1))
  # dual_temp = solve.QP(D, dvec, Amat, bvec, meq = n_class - 1)
   ###
  alpha = dual$solution
  alpha[alpha < 0] = 0

  alpha_mat = matrix(alpha, nrow = n_l, ncol = n_class)
  # alpha_mat[y_index][alpha_mat[y_index] > gamma] = gamma

  # for (j in 1:n_class) {
  #   alpha_mat[y != j, j][alpha_mat[y != j, j] > (1 - gamma)] = (1 - gamma)
  # }

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

  # alpha_vec = as.vector(alpha_mat)

  # cmat = matrix(0, n, n_class - 1)
  # for (k in 1:(n_class - 1)) {
  #   cmat[, k] = inv_KLK %*% K %*% t(J) %*% t(Hmatj[[k]]) %*% alpha_vec
  # }

  cmat = matrix(0, n, n_class - 1)
  for (k in 1:(n_class - 1)) {
    cmat[, k] = inv_K_KLK %*% t(Hmatj[[k]]) %*% alpha
    # cmat[, k] = inv_K_KLK %*% t(JK) %*% t(Hmatj[[k]]) %*% alpha
  }

  # find b vector using LP
  Kcmat = (JK %*% cmat) %*% W

  # table(y, apply(Kcmat, 1, which.max))


  alp_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  alp_temp[y_index] = gamma

  alp = c(as.vector(alp_temp), rep(0, 2 * (n_class - 1)))

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
  # compute the fitted values
  fit = (matrix(W_c0vec, nrow = n_l, ncol = n_class, byrow = T) + Kcmat)
  fit_class = levs[apply(fit, 1, which.max)]
  if (attr(levs, "type") == "factor") {fit_class = factor(fit_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {fit_class = as.numeric(fit_class)}
  if (attr(levs, "type") == "integer") {fit_class = as.integer(fit_class)}
  # table(y, fit_class)

  # Return the output
  out$alpha = alpha_mat
  out$beta = cmat
  out$beta0 = c0vec
  out$fit = fit
  out$fit_class = fit_class
  out$n_l = n_l
  out$n_u = n_u
  out$n_class = n_class
  out$levels = levs
  return(out)
}
