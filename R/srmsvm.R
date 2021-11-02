srmsvm = function(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                 kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                 scale = TRUE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep.srmsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                          lambda_seq = lambda_seq, theta = NULL,
                          kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  cat("Fit theta-step \n")
  thetastep_fit = thetastep.srmsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.srmsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                              lambda_seq = lambda_seq, theta = thetastep_fit$opt_theta,
                              kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)

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
  out$n_class = opt_cstep_fit$n_class
  class(out) = "srmsvm"
  return(out)
}

predict.srmsvm = function(object, newx = NULL, newK = NULL)
{
  model = object$opt_model
  cmat = model$cmat
  c0vec = model$c0vec
  levs = model$levels

  # if (object$scale) {
  #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  # }

  if (is.null(newK)) {
    new_anova_K = make_anovaKernel(newx, object$x, kernel = object$kernel, kparam = object$kparam)
    newK = combine_kernel(new_anova_K, theta = object$opt_theta)
    # newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = model$n_class, byrow = T) + (newK %*% cmat))
  pred_class = apply(pred_y, 1, which.max)

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}


cstep.srmsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                       lambda_seq = 2^{seq(-10, 10, length.out = 100)}, theta = NULL,
                       kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                       scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  p = ncol(x)

  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)

  lambda_seq = sort(lambda_seq, decreasing = FALSE)

  n = NROW(x)
  center = rep(0, p)
  scaled = rep(1, p)

  if (scale) {
    x = scale(x)
    center = attr(x, "scaled:center")
    scaled = attr(x, "scaled:scale")
  }

  anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)
  if (is.null(theta)) {
    theta = rep(1, anova_K$numK)
  }
  K = combine_kernel(anova_K, theta)

  # Combination of hyper-parameters

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL
    ran = NULL

    valid_anova_K = make_anovaKernel(valid_x, x, kernel = kernel, kparam = kparam)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:length(lambda_seq),
                        function(j) {
                          error = try({
                            msvm_fit = rmsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda_seq[j], ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.rmsvm_compact(msvm_fit, newK = valid_K)
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
    valid_err = sapply(fold_err, "[[", "error")
    # model_list[[1]] = lapply(fold_err, "[[", "fit_model")

    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
  } else {
    ran = data_split(y, nfolds)
    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_seq), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i_cv in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      omit = ran == i_cv
      x_train = x[!omit, ]
      y_train = y[!omit]
      x_valid = x[omit, ]
      y_valid = y[omit]

      subanova_K = make_anovaKernel(x_train, x_train, kernel, kparam)
      subK = combine_kernel(subanova_K, theta)

      subanova_K_valid = make_anovaKernel(x_valid, x_train, kernel, kparam)
      subK_valid = combine_kernel(subanova_K_valid, theta)

      fold_err = mclapply(1:length(lambda_seq),
                          function(j) {
                            error = try({
                              msvm_fit = rmsvm_compact(K = subK, y = y_train, gamma = gamma, lambda = lambda_seq[j], ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.rmsvm_compact(msvm_fit, subK_valid)
                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {
                                # 수정 필요 y_valid가 factor나 character일 경우
                                err = ramsvm_hinge(y_valid, pred_val$pred_value, k = k, gamma = gamma)
                              }
                            } else {
                              msvm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = msvm_fit))
                          }, mc.cores = nCores)

      valid_err_mat[i_cv, ] = sapply(fold_err, "[[", "error")
    }
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
  }
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$gamma = gamma
  out$theta = theta
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = kparam
  out$scale = scale
  out$criterion = criterion
  out$nfolds = nfolds
  out$fold_list = ran
  if (optModel) {
    # anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)
    # K = combine_kernel(anova_K, theta)
    opt_model = rmsvm_compact(K = K, y = y, gamma = gamma, lambda = out$opt_param["lambda"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "srmsvm"
  return(out)
}

thetastep.srmsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$kparam
  x = object$x
  y = object$y
  gamma = object$gamma
  theta = object$theta
  nfolds = object$nfolds
  fold_list = object$fold_list

  valid_x = object$valid_x
  valid_y = object$valid_y

  anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)

  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, theta)
    opt_model = rmsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  } else {
    opt_model = object$opt_model
  }

  if (!is.null(valid_x) & !is.null(valid_y)) {
    valid_anova_K = make_anovaKernel(valid_x, x, kernel = kernel, kparam = kparam)

    init_model = opt_model

    fold_err = mclapply(1:length(lambda_theta_seq),
                        function(j) {
                          error = try({
                            theta = findtheta.srmsvm(y = y, anova_kernel = anova_K, gamma = gamma,
                                                      cmat = init_model$cmat, c0vec = init_model$c0vec,
                                                      lambda = lambda, lambda_theta = lambda_theta_seq[j], ...)
                            if (isCombined) {
                              subK = combine_kernel(anova_K, theta)
                              init_model = rmsvm_compact(K = subK, y = y, gamma = gamma, lambda = lambda, ...)
                            }
                          })

                          if (!inherits(error, "try-error")) {
                            valid_subK = combine_kernel(valid_anova_K, theta)
                            pred_val = predict.rmsvm_compact(init_model, newK = valid_subK)

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
    valid_err = sapply(fold_err, "[[", "error")
    theta_seq = sapply(fold_err, "[[", "theta")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_theta = theta_seq[, opt_ind]
    opt_valid_err = min(valid_err)
  } else {

    valid_err_mat = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i_cv in 1:nfolds) {
      cat(nfolds, "- fold CV :", i / nfolds * 100, "%", "\r")
      omit = fold_list == i_cv
      x_train = x[!omit, ]
      y_train = y[!omit]
      x_valid = x[omit, ]
      y_valid = y[omit]

      subanova_K = make_anovaKernel(x_train, x_train, kernel, kparam)
      subK = combine_kernel(subanova_K, theta)
      subanova_K_valid = make_anovaKernel(x_valid, x_train, kernel, kparam)

      init_model = rmsvm_compact(K = subK, y = y_train, gamma = gamma, lambda = lambda, ...)
      cmat = init_model$cmat
      c0vec = init_model$c0vec

      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = findtheta.srmsvm(y = y_train, anova_kernel = subanova_K, gamma = gamma, cmat = cmat, c0vec = c0vec,
                                                        lambda = lambda, lambda_theta = lambda_theta_seq[j])
                              if (isCombined) {
                                subK = combine_kernel(subanova_K, theta)
                                init_model = rmsvm_compact(K = subK, y = y_train, gamma = gamma, lambda = lambda, ...)
                              }
                            })

                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.rmsvm_compact(init_model, newK = subK_valid)

                              if (criterion == "0-1") {
                                acc = sum(y_valid == pred_val$class) / length(y_valid)
                                err = 1 - acc
                              } else {
                                # 수정 필요 y_valid가 factor나 character일 경우
                                err = rmsvm_hinge(y_valid, pred_val$pred_value, k = k, gamma = gamma)
                              }
                            } else {
                              err = Inf
                              theta = rep(0, anova_K$numK)
                            }
                            return(list(error = err, theta = theta))
                          }, mc.cores = nCores)
      valid_err_mat[i_cv, ] = sapply(fold_err, "[[", "error")
    }
    valid_err = colMeans(valid_err_mat)
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(valid_err)

    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = findtheta.srmsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
                                                            lambda = lambda, lambda_theta = lambda_theta_seq[j])
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

  if (optModel) {
    optK = combine_kernel(anova_K, opt_theta)
    opt_model = rmsvm_compact(K = optK, y = y, gamma = gamma, lambda = lambda, ...)
    out$opt_model = opt_model
  }
  class(out) = "srmsvm"
  return(out)
}


findtheta.srmsvm = function(y, anova_kernel, gamma, cmat, c0vec, lambda, lambda_theta, eig_tol_D = 0, epsilon_D = 1e-8)
{

  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }

  if (lambda_theta <= 0) {
    theta = rep(1, anova_kernel$numK)
    return(theta)
  }

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)
  n_class = length(levs)

  # standard LP form :
  # min a^T x , subject to A1x <= a1
  n = length(y_int)
  y_index = cbind(1:n, y_int)

  a_tmp = matrix(gamma / n, nrow = n, ncol = n_class)
  a_tmp[y_index] = (1 - gamma) / n
  a = matrix(a_tmp, ncol = 1)

  # initialize M
  M = matrix(rep(0, anova_kernel$numK), ncol = 1)
  # calculate M
  for (d in 1:anova_kernel$numK) {
    for (j in 1:n_class) {
      M[d] = (M[d] + t(cmat[, j]) %*% anova_kernel$K[[d]] %*% cmat[, j])
    }
    M[d] = (lambda / 2 * M[d] + (lambda_theta))
  }
  a = rbind(a, M)
  # calculate N matrix
  for (d in 1:anova_kernel$numK) {
    K = anova_kernel$K[[d]]
    for (j in 1:n_class) {
      if (j == 1) {
        temp_N = K %*% cmat[, j]
      } else {
        temp_N = rbind(temp_N, K %*% cmat[, j])
      }
    }
    if (d == 1) {
      N = temp_N
    } else {
      N = cbind(N, temp_N)
    }
  }

  # constraints
  sign_mat = matrix(1, n, n_class)
  sign_mat[y_index] = -1

  # A matrix
  I_nonzeroIndex = diag(1, n * n_class)
  N_nonzeroIndex = as.matrix(N) * as.vector(sign_mat)
  A_theta = cbind(matrix(0, anova_kernel$numK, n * n_class), diag(-1, anova_kernel$numK))
  A_ineq = rbind(cbind(I_nonzeroIndex, -N_nonzeroIndex), A_theta)

  bb = matrix(rep(c0vec, n), n, n_class, byrow = TRUE)
  bb_j = 1 + bb
  bb_yi = (n_class - 1) - bb
  bb_j[y_index] = bb_yi[y_index]

  b_nonzeroIndex = as.vector(bb_j)
  b_ineq = rbind(as.matrix(b_nonzeroIndex), matrix(-1, anova_kernel$numK, 1))

  # constraint directions
  const_dir = matrix(rep(">=", nrow(matrix(b_ineq))))
  # find solution by LP
  lp = lp("min", objective.in = a, const.mat = A_ineq, const.dir = const_dir,
          const.rhs = b_ineq)$solution
  # find the theta vector only from the solution
  theta = cbind(matrix(0, anova_kernel$numK, n * n_class),
                diag(1, anova_kernel$numK)) %*% matrix(lp, ncol = 1)
  return(as.vector(theta))
}

































