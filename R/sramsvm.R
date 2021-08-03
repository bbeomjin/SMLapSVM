# sramsvm = function(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
#                   lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
#                   kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
#                   scale = TRUE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, verbose = 1, ...)
# {
#   out = list()
#   call = match.call()
#   kernel = match.arg(kernel)
#   criterion = match.arg(criterion)
#
#   cat("Fit c-step \n")
#
#   cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
#                            lambda_seq = lambda_seq, theta = NULL, kernel = kernel, kparam = kparam,
#                            criterion = criterion, optModel = FALSE, nCores = nCores, ...)
#
#   cat("Fit theta-step \n")
#
#   if (optModel) {
#     thetastep_opt = FALSE
#   } else {
#     thetastep_opt = TRUE
#   }
#
#   thetastep_fit = thetastep.sramsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined,
#                                     optModel = thetastep_opt, nCores = nCores, ...)
#
#   if (verbose == 1)
#
#   cat("Fit c-step \n")
#   opt_cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
#                                lambda_seq = lambda_seq, theta = thetastep_fit$opt_theta,
#                                kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
#
#   out$opt_param = opt_cstep_fit$opt_param
#   out$opt_valid_err = opt_cstep_fit$opt_valid_err
#   out$cstep_valid_err = opt_cstep_fit$valid_err
#   out$theta_valid_err = thetastep_fit$valid_err
#   out$opt_model = opt_cstep_fit$opt_model
#   out$kernel = kernel
#   out$kparam = opt_cstep_fit$opt_param["kparam"]
#   out$opt_theta = thetastep_fit$opt_theta
#   out$theta = thetastep_fit$theta
#   out$x = x
#   out$y = y
#   out$n_class = opt_cstep_fit$n_class
#   class(out) = "sramsvm"
#   return(out)
# }

# predict.sramsvm = function(object, newx = NULL, newK = NULL)
# {
#   model = object$opt_model
#   cmat = model$cmat
#   c0vec = model$c0vec
#   levs = model$levels
#
#   W = XI_gen(model$n_class)
#
#   # if (object$scale) {
#   #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
#   # }
#
#   if (is.null(newK)) {
#     new_anova_K = make_anovaKernel(newx, object$x, kernel = object$kernel, kparam = object$kparam)
#     newK = combine_kernel(new_anova_K, theta = object$opt_theta)
#     # newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
#     # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
#   }
#
#   W_c0 = drop(t(c0vec) %*% W)
#
#   pred_y = matrix(W_c0, nrow = nrow(newK), ncol = model$n_class, byrow = T) + ((newK %*% cmat) %*% W)
#   pred_class = levs[apply(pred_y, 1, which.max)]
#
#   if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
#   if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
#   if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}
#
#   return(list(class = pred_class, pred_value = pred_y))
# }


cstep.sramsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                        lambda_seq = 2^{seq(-10, 10, length.out = 100)}, theta = NULL,
                        kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                        criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  n = NROW(x)
  p = ncol(x)

  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)

  lambda_seq = sort(lambda_seq, decreasing = FALSE)

  # Combination of hyper-parameters

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)

    if (is.null(theta)) {
      theta = rep(1, anova_K$numK)
    }

    K = combine_kernel(anova_K, theta)

    valid_anova_K = make_anovaKernel(valid_x, x, kernel = kernel, kparam = kparam)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:length(lambda_seq),
                        function(j) {
                          error = try({
                            msvm_fit = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda_seq[j], ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.ramsvm_compact(msvm_fit, newK = valid_K)$class
                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val) / length(valid_y)
                              err = 1 - acc
                            } else {
                              # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
                            }
                          } else {
                            smsvm_fit = NULL
                            err = Inf
                          }
                          return(list(error = err, fit_model = msvm_fit))
                        }, mc.cores = nCores)
    valid_err = sapply(fold_err, "[[", "error")
    # model_list[[1]] = lapply(fold_err, "[[", "fit_model")

    opt_ind = max(valid_err == min(valid_err))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(valid_err)
  } else {
    ran = data_split(y, nfolds)
    valid_err = matrix(NA, nrow = nfolds, ncol = length(lambda_seq), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]

      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      if (is.null(theta)) {
        theta = rep(1, subanova_K$numK)
      }
      subK = combine_kernel(subanova_K, theta)

      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)
      subK_valid = combine_kernel(subanova_K_valid, theta)

      fold_err = mclapply(1:length(lambda_seq),
                          function(j) {
                            error = try({
                              msvm_fit = ramsvm_compact(K = subK, y = train_y, gamma = gamma, lambda = lambda_seq[j], ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.ramsvm_core(msvm_fit, subK_valid)
                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                # 수정 필요 valid_y가 factor나 character일 경우
                                err = ramsvm_hinge(valid_y, pred_val$pred_value, k = k, gamma = gamma)
                              }
                            } else {
                              msvm_fit = NULL
                              err = Inf
                            }
                            return(list(error = err, fit_model = msvm_fit))
                          }, mc.cores = nCores)

      valid_err[i_cv, ] = sapply(fold_err, "[[", "error")
    }
    mean_valid_err = colMeans(valid_err)
    opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
    opt_param = c(lambda = lambda_seq[opt_ind])
    opt_valid_err = min(mean_valid_err)
    valid_x = NULL
    valid_y = NULL
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
  out$criterion = criterion
  out$nfolds = nfolds
  if (optModel) {
    anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)
    K = combine_kernel(anova_K, theta)
    opt_model = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = opt_param["lambda"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "sramsvm"
  return(out)
}

thetastep.sramsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, optModel = FALSE, nCores = 1, ...)
{
  out = list()
  call = match.call()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  kernel = object$kernel
  kparam = object$kparam
  criterion = object$criterion
  n_class = object$n_class
  x = object$x
  y = object$y
  gamma = object$gamma
  theta = object$theta

  valid_x = object$valid_x
  valid_y = object$valid_y

  anova_K = make_anovaKernel(x, x, kernel = kernel, kparam = kparam)
  valid_anova_K = make_anovaKernel(valid_x, x, kernel = kernel, kparam = kparam)

  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, theta)
    opt_model = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  } else {
    opt_model = object$opt_model
  }

  if (!is.null(valid_x) & !is.null(valid_y)) {

    init_model = opt_model
    cmat = init_model$cmat
    c0vec = init_model$c0vec

    fold_err = mclapply(1:length(lambda_theta_seq),
                        function(j) {
                          error = try({
                            theta = findtheta.sramsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = cmat, c0vec = c0vec,
                                                      lambda = lambda, lambda_theta = lambda_theta_seq[j], ...)
                            if (isCombined) {
                              subK = combine_kernel(anova_K, theta)
                              init_model = ramsvm_compact(K = subK, y = y, gamma = gamma, lambda = lambda, ...)
                            }
                          })

                          if (!inherits(error, "try-error")) {
                            valid_subK = combine_kernel(valid_anova_K, theta)
                            pred_val = predict.ramsvm_compact(init_model, newK = valid_subK)$class

                            if (criterion == "0-1") {
                              acc = sum(valid_y == pred_val) / length(valid_y)
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

    ran = data_split(y, nfolds)
    valid_err = matrix(NA, nrow = nfolds, ncol = length(lambda_theta_seq), dimnames = list(paste0("Fold", 1:nfolds)))

    for (i_cv in 1:nfolds) {
      omit = ran == i_cv
      train_x = x[!omit, ]
      train_y = y[!omit]
      valid_x = x[omit, ]
      valid_y = y[omit]

      subanova_K = make_anovaKernel(train_x, train_x, kernel, kparam)
      subK = combine_kernel(subanova_K, theta)
      subanova_K_valid = make_anovaKernel(valid_x, train_x, kernel, kparam)

      init_model = ramsvm_compact(K = subK, y = train_y, gamma = gamma, lambda = lambda, ...)
      cmat = init_model$cmat
      c0vec = init_model$c0vec

      fold_err = mclapply(1:length(lambda_theta_seq),
                          function(j) {
                            error = try({
                              theta = findtheta.sramsvm(y = train_y, anova_kernel = subanova_K, gamma = gamma, cmat = cmat, c0vec = c0vec,
                                                        lambda = lambda, lambda_theta = lambda_theta_seq[j])
                              if (isCombined) {
                                subK = combine_kernel(subanova_K, theta)
                                init_model = ramsvm_compact(K = subK, y = train_y, gamma = gamma, lambda = lambda, ...)
                              }
                            })

                            if (!inherits(error, "try-error")) {
                              subK_valid = combine_kernel(subanova_K_valid, theta)
                              pred_val = predict.ramsvm_core(init_model, newK = subK_valid)

                              if (criterion == "0-1") {
                                acc = sum(valid_y == pred_val$class) / length(valid_y)
                                err = 1 - acc
                              } else {
                                # 수정 필요 valid_y가 factor나 character일 경우
                                err = ramsvm_hinge(valid_y, pred_val$pred_value, k = k, gamma = gamma)
                              }
                            } else {
                              err = Inf
                              theta = rep(0, anova_K$numK)
                            }
                            return(list(error = err, theta = theta))
                          }, mc.cores = nCores)
      valid_err[i_cv, ] = sapply(fold_err, "[[", "error")
    }
    mean_valid_err = colMeans(valid_err)
    opt_ind = max(which(mean_valid_err == min(mean_valid_err)))
    opt_lambda_theta = lambda_theta_seq[opt_ind]
    opt_valid_err = min(mean_valid_err)

    theta_seq_list = mclapply(1:length(lambda_theta_seq),
                              function(j) {
                                error = try({
                                  theta = findtheta.sramsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = opt_model$cmat, c0vec = opt_model$c0vec,
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
    opt_model = ramsvm_compact(K = optK, y = y, gamma = gamma, lambda = lambda, ...)
    out$opt_model = opt_model
  }
  class(out) = "sramsvm"
  return(out)
}


findtheta.sramsvm = function(y, anova_kernel, gamma, cmat, c0vec, lambda, lambda_theta, eig_tol_D = 0, epsilon_D = 1e-8)
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

  n = length(y_int)
  y_index = cbind(1:n, y_int)

  # c0vec = as.matrix(c0vec)

  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, nobs = n, y = y_int)

  # calculate the 'a' matrix
  a_tmp = matrix(gamma / n, nrow = n, ncol = n_class)
  a_tmp[y_index] = (1 - gamma) / n
  a = matrix(a_tmp, ncol = 1)

  # initialize M
  M = matrix(rep(0, anova_kernel$numK), ncol = 1)

  # calculate M
  for (d in 1:anova_kernel$numK) {
    for (j in 1:(n_class - 1)) {
      M[d] = (M[d] + t(cmat[, j]) %*% anova_kernel$K[[d]] %*% cmat[, j])
    }
    M[d] = (lambda / 2 * M[d] + (lambda_theta))
  }
  a = rbind(a, M)

  Y_code = XI_gen(n_class)

  # calculate N matrix
  for (d in 1:anova_kernel$numK) {
    K = anova_kernel$K[[d]]
    for (j in 1:n_class) {
      if(j == 1) {
        temp_N = matrix(rowSums(K %*% cmat * matrix(Y_code[, j], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE)))
      } else {
        temp_N = rbind(temp_N, matrix(rowSums(K %*% cmat * matrix(Y_code[, j], nrow = nrow(K), ncol = ncol(cmat), byrow = TRUE))))
      }
    }
    if (d == 1) {
      N = temp_N
    } else {
      N = cbind(N, temp_N)
    }
  }

  sign_mat = matrix(1, n, n_class)
  sign_mat[y_index] = -1

  # constraints

  # A matrix
  I_nonzeroIndex = diag(1, n * n_class)
  N_nonzeroIndex = as.matrix(N) * as.vector(sign_mat)
  A_theta = cbind(matrix(0, anova_kernel$numK, n * n_class),
                  diag(-1, anova_kernel$numK))
  A_ineq = rbind(cbind(I_nonzeroIndex, -N_nonzeroIndex), A_theta)

  bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * c0vec[x]))
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(rowSums(t(Y_code) * matrix(c0vec, nrow = ncol(Y_code), ncol = nrow(Y_code), byrow = TRUE)), nrow = n, ncol = n_class, byrow = TRUE)
  bb_j[y_index] = bb_yi

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

































