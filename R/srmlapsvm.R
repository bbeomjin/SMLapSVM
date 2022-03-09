srmlapsvm = function(x = NULL, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                    lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                    adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                    kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                    scale = FALSE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, verbose = 0, ...)
{
  out = list()

  fun_ind = ifelse(gamma == 0.0, 1, 2)

  cstep_fun = switch(fun_ind,
                     cstep.smlapsvm,
                     cstep.srmlapsvm)

  thetastep_fun = switch(fun_ind,
                          thetastep.smlapsvm,
                          thetastep.srmlapsvm)

  cat("Fit c-step \n")
  # cstep_fit = cstep.srmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
  #                   lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
  #                   adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
  #                   kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  cstep_args = list(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq,
                    theta = NULL, fold_theta = NULL,
                    adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                    kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  if (fun_ind == 2) {cstep_args$gamma = gamma}
  cstep_fit = do.call(cstep_fun, cstep_args)

  cat("Fit theta-step \n")
  thetastep_fit = thetastep_fun(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)


  cat("Fit c-step \n")
  opt_cstep_args = list(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                        lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq,
                        theta = thetastep_fit$opt_theta, fold_theta = thetastep_fit$opt_fold_theta,
                        adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                        kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
  # opt_cstep_fit = cstep.srmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
  #                       lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = thetastep_fit$opt_theta,
  #                       adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
  #                       kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
  if (fun_ind == 2) {opt_cstep_args$gamma = gamma}
  opt_cstep_fit = do.call(cstep_fun, opt_cstep_args)

  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
	cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
	cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
  }

  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$cstep_valid_err = opt_cstep_fit$valid_err
  out$theta_valid_err = thetastep_fit$valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param["kparam"]
  out$opt_theta = thetastep_fit$opt_theta
  out$theta = thetastep_fit$theta
  out$x = x
  out$y = y
  out$ux = ux
  out$gamma = gamma
  out$n_class = opt_cstep_fit$n_class
  class(out) = "srmlapsvm"
  return(out)
}

predict.srmlapsvm = function(object, newx = NULL, newK = NULL)
{
  model = object$opt_model
  cmat = model$cmat
  c0vec = model$c0vec
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

  pred_y = (matrix(c0vec, nrow = nrow(newK), ncol = model$n_class, byrow = T) + (newK %*% cmat))
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}


cstep.srmlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                 theta = NULL, fold_theta = NULL,
                 adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                 kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                 scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  p = ncol(x)

  lambda_seq = as.numeric(lambda_seq)
  lambda_I_seq = as.numeric(lambda_I_seq)
  kparam = as.numeric(kparam)

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  lambda_I_seq = sort(lambda_I_seq, decreasing = FALSE)
  # kparam = sort(kparam, decreasing = FALSE)


  # 추후, 커널에 맞게 theta의 길이 조절
  if (is.null(theta)) {
    theta = rep(1, p)
  }

  if (is.null(fold_theta)) {
    fold_theta = rep(list(rep(1, p)), nfolds)
  }

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

    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = kparam)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            msvm_fit = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                         lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.rmlapsvm_compact(msvm_fit, newK = valid_K)
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
                              msvm_fit = srmlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta_train, y = y_train, gamma = gamma,
                                                            lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.rmlapsvm_compact(msvm_fit, newK = subK_valid)

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
    # valid_err = round(colMeans(valid_err_mat), 8)
    valid_err = colMeans(valid_err_mat)
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
    opt_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                  lambda = out$opt_param["lambda"], lambda_I = out$opt_param["lambda_I"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  return(out)
}


thetastep.srmlapsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
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
    opt_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = init_theta, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
  } else {
    opt_model = object$opt_model
  }

  if (!is.null(valid_x) & !is.null(valid_y)) {

    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = kparam)

    init_model = opt_model

    fold_err = mclapply(1:length(lambda_theta_seq),
                        function(j) {
                          error = try({
                            theta = find_theta.srmlapsvm(y = y, gamma = gamma, anova_kernel = anova_K, L = L,
                                                         cmat = init_model$cmat, c0vec = init_model$c0vec,
                                                         lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
                            if (isCombined) {
                              # subK = combine_kernel(anova_K, theta)
                              init_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                             lambda = lambda, lambda_I = lambda_I, ...)
                            }
                          })

                          if (!inherits(error, "try-error")) {
                            valid_subK = combine_kernel(valid_anova_K, theta)
                            pred_val = predict.rmlapsvm_compact(init_model, newK = valid_subK)

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
    # valid_err = round(sapply(fold_err, "[[", "error"), 8)
    valid_err = sapply(fold_err, "[[", "error")
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
    #                               theta = find_theta.srmlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
    #                                                            gamma = gamma, lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
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

      init_model = srmlapsvm_compact(anova_K = subanova_K, L = L_train, theta = init_theta_train, y = y_train,
                                      lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
      cmat = init_model$cmat
      c0vec = init_model$c0vec

      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = find_theta.srmlapsvm(y = y_train, anova_kernel = subanova_K, L = L_train, cmat = cmat, c0vec = c0vec,
                                                            gamma = gamma, lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)

                              if (isCombined) {
                                init_model = srmlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta, y = y_train,
                                                                lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
                              }
                            })

                            # error = try({
                            #   theta = theta_seq_list[[j]]
                            #   if (isCombined) {
                            #     init_model = srmlapsvm_compact(anova_K = subanova_K, L = L_train, theta = theta, y = y_train,
                            #                                    lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
                            #   }
                            # })

                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.rmlapsvm_compact(init_model, newK = subK_valid)

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
    # valid_err = round(colMeans(valid_err_mat), 8)
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)

    opt_fold_theta = lapply(fold_theta, function(x) return(x[, opt_ind]))

    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = find_theta.srmlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
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
    opt_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = opt_theta, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
    out$opt_model = opt_model
  }
  class(out) = "srmlapsvm"
  return(out)
}


find_theta.srmlapsvm = function(y, gamma, anova_kernel, L, cmat, c0vec, lambda, lambda_I, lambda_theta = 1,
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

  # n = NROW(cmat)
  n_l = length(y_int)
  n_u = n - n_l

  # Y = class_code(y, k = n_class)
  y_index = cbind(1:n_l, y_int)

  Dmat = numeric(anova_kernel$numK)
  dvec = numeric(anova_kernel$numK)
  A_mat = NULL

  for(j in 1:anova_kernel$numK) {
    temp_D = 0
    temp_d = 0
    temp_A = NULL
    for (q in 1:ncol(cmat)) {
      cvec = cmat[, q]
      # KLK_temp = anova_kernel_orig$K[[j]] %*% L %*% anova_kernel_orig$K[[j]]
      # diag(KLK_temp) = diag(KLK_temp) + max(abs(KLK_temp)) * epsilon_I

      temp_D = temp_D + n_l * lambda_I / (1 * n^2) * t(cvec) %*% anova_kernel$K[[j]] %*% L %*% anova_kernel$K[[j]] %*% cvec
      temp_d = temp_d + n_l * lambda / 2 * t(cvec) %*% anova_kernel$K[[j]] %*% cvec + n_l * lambda_theta
      temp_A = rbind(temp_A, (anova_kernel$K[[j]][1:n_l, ] %*% cvec))
    }
    Dmat[j] = temp_D
    dvec[j] = temp_d
    A_mat = cbind(A_mat, temp_A)
  }

  max_D = max(abs(Dmat))
  Dmat = c(Dmat, c(rep(0, n_l * n_class)))
  Dmat = diag(Dmat)
  # Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)
  # diag(Dmat) = diag(Dmat) + max_D * epsilon_D
  diag(Dmat) = diag(Dmat) + max_D * epsilon_D

  # Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)

  # Dmat = Dmat / max_D

  dvec_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  dvec_temp[y_index] = gamma

  # dvec_temp = as.vector(Y)
  # dvec_temp[dvec_temp == 1] = 0
  # dvec_temp[dvec_temp < 0] = 1
  dvec = c(dvec, as.vector(dvec_temp))
  # dvec = dvec
  # dvec = dvec / max_D

  # solve QP

  # diag(Dmat) = diag(Dmat) + epsilon_D

  m_index = matrix(1:(n_l * n_class), ncol = n_class)[y_index]
  A_mat[m_index, ] = -A_mat[m_index, ]
  A_mat = cbind(-A_mat, diag(1, n_l * n_class))
  A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)

  bb = c0vec[y_int]
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(c0vec, nrow = n_l, ncol = n_class, byrow = TRUE)
  bb_j[y_index] = bb_yi
  bvec = c(as.vector(bb_j), rep(0, anova_kernel$numK + n_l * n_class), rep(-1, anova_kernel$numK))

  theta_sol = solve.QP(Dmat, -dvec, t(A_mat), bvec, meq = 0, factorized = FALSE)$solution
  theta = theta_sol[1:anova_kernel$numK]
  theta[theta < 1e-6] = 0
  theta = round(theta, 6)

  return(theta)
}



srmlapsvm_compact = function(anova_K, L, theta, y, gamma = 0.5, lambda, lambda_I,
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

  # anova_K_orig = anova_K
  # anova_K$K = lapply(anova_K$K, function(x) {
  #   diag(x) = diag(x) + max(abs(diag(x))) * epsilon_I
  #   return(x)
  # })

  K = combine_kernel(anova_K, theta = theta)
  n = nrow(K)
  # K = (K + t(K)) / 2

  if (sum(K) == 0) {
    diag(K) = 1
  }

  # n = nrow(K)
  n_l = length(y_int)
  n_u = n - n_l
  qp_dim = n_l * n_class

  code_mat = code_rmsvm(y_int)
  In = code_mat$In
  vmatj = code_mat$vmatj
  umatj = code_mat$umatj
  Hmatj = code_mat$Hmatj
  y_index = code_mat$y_index

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

  # K = fixit(K, epsilon = eig_tol_I)
  # KLK = fixit(KLK, epsilon = eig_tol_I)

  lambda_K = n_l * lambda * K
  lambda_K = fixit(lambda_K, epsilon = eig_tol_I)
  # diag(lambda_K) = diag(lambda_K) + epsilon_I
  # diag(lambda_K) = diag(lambda_K) + max(abs(diag(lambda_K))) * epsilon_I

  lambda_KLK = n_l * lambda_I / n^2 * KLK
  lambda_KLK = fixit(lambda_KLK, epsilon = eig_tol_I)
  # lambda_KLK = fixit(lambda_KLK)
  # diag(lambda_KLK) = diag(lambda_KLK) + epsilon_I
  # diag(lambda_KLK) = diag(lambda_KLK) + max(abs(diag(lambda_KLK))) * epsilon_I

  K_KLK = lambda_K + lambda_KLK
  # K_KLK = n_l * lambda * K + n_l * lambda_I / n^2 * KLK
  # K_KLK = (K_KLK + t(K_KLK)) / 2
  # K_KLK = fixit(K_KLK, epsilon = eig_tol_I)
  # K_KLK = fixit(K_KLK)
  # diag(K_KLK) = diag(K_KLK) + max(abs(diag(K_KLK))) * epsilon_I
  # diag(K_KLK) = diag(K_KLK) + max(abs(diag(K_KLK))) * epsilon_I

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
  # Q = fixit(Q, epsilon = eig_tol_D)
  # diag(Q) = diag(Q) + epsilon_D

  # Compute Q = K x inv_LK
  D = 0
  Amat = matrix(0, (2 * qp_dim + n_class), qp_dim)
  for (j in 1:n_class) {
    D = D + t(Hmatj[[j]]) %*% Q %*% Hmatj[[j]]
    Amat[j, ] = rep(1, n_l) %*% Hmatj[[j]]
  }
  D = (D + t(D)) / 2
  D = fixit(D, epsilon = eig_tol_D)
  # D = fixit(D)
  # D = fixit2(D, epsilon = 0)
  max_D = max(abs(diag(D)))
  # D = D / max_D
  diag(D) = diag(D) + max_D * epsilon_D
  # diag(D) = diag(D) + epsilon_D

  # D = nearPD(D, eig.tol = rel_eig_tol)$mat
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

  cmat = matrix(0, n, n_class)
  for (k in 1:n_class) {
    cmat[, k] = inv_K_KLK %*% Hmatj[[k]] %*% alpha
    # cmat[, k] = inv_K_KLK %*% t(JK) %*% Hmatj[[k]] %*% alpha
  }

  # find b vector using LP
  Kcmat = JK %*% cmat

  alp_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  alp_temp[y_index] = gamma

  alp = c(as.vector(alp_temp), rep(0, 2 * n_class))

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

  ############################################################################

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
















































