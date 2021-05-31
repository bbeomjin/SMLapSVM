sramsvm = function(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                  lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                  kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                  scale = TRUE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                           lambda_seq = lambda_seq, theta = NULL,
                           kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  cat("Fit theta-step \n")
  theta_step_fit = theta_step.sramsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.sramsvm(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                               lambda_seq = lambda_seq, theta = theta_step_fit$opt_theta,
                               kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)

  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$cstep_valid_err = opt_cstep_fit$valid_err
  out$theta_valid_err = theta_step_fit$valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param["kparam"]
  out$opt_theta = theta_step_fit$opt_theta
  out$theta = theta_step_fit$theta
  out$x = x
  out$y = y
  out$n_class = opt_cstep_fit$n_class
  class(out) = "sramsvm"
  return(out)
}

predict.sramsvm = function(object, newx = NULL, newK = NULL)
{
  model = object$opt_model
  kernel = object$kernel
  kparam = object$kparam
  beta = model$beta
  beta0 = model$beta0
  n_class = object$n_class
  W = XI_gen(n_class)

  # if (object$scale) {
  #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  # }

  if (is.null(newK)) {
    new_anova_K = make_anovaKernel(newx, object$x, kernel = list(type = object$kernel, par = object$kparam))
    newK = combine_kernel(new_anova_K, theta = object$opt_theta)
    # newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  W_beta0 = drop(t(beta0) %*% W)

  fit = matrix(W_beta0, nrow = nrow(newK), ncol = n_class, byrow = T) + ((newK %*% beta) %*% W)
  pred_y = apply(fit, 1, which.max)
  return(list(class = pred_y, pred_value = fit))
}


cstep.sramsvm = function(x, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                        lambda_seq = 2^{seq(-10, 10, length.out = 100)}, theta = NULL,
                        kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                        scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  p = ncol(x)

  lambda_seq = as.numeric(lambda_seq)
  kparam = as.numeric(kparam)

  if (is.null(theta)) {
    theta = rep(1, p)
  }

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  kparam = sort(kparam, decreasing = FALSE)

  # Combination of hyper-parameters

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    n = NROW(x)
    # The number of classes
    n_class = length(unique(y))

    center = rep(0, p)
    scaled = rep(1, p)

    if (scale) {
      x = scale(x)
      center = attr(x, "scaled:center")
      scaled = attr(x, "scaled:scale")
    }

    valid_err_mat = matrix(NA, nrow = length(kparam), ncol = length(lambda_seq))

    for (i in 1:length(kparam)) {
      par = kparam[i]

      kernel_list = list(type = kernel, par = par)
      anova_K = make_anovaKernel(x, x, kernel = kernel_list)
      K = combine_kernel(anova_K, theta)

      valid_anova_K = make_anovaKernel(valid_x, x, kernel = kernel_list)
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
      valid_err_mat[i, ] = valid_err
    }
    opt_ind = which(valid_err_mat == min(valid_err_mat), arr.ind = TRUE)
    opt_ind = opt_ind[order(opt_ind[, 1], opt_ind[, 2], decreasing = c(FALSE, TRUE))[1], ]
    opt_param = c(lambda = lambda_seq[opt_ind[2]], kparam = kparam[opt_ind[1]])
    opt_valid_err = min(valid_err_mat)
  }
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$gamma = gamma
  out$theta = theta
  out$n_class = n_class
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$kparam = opt_param["kparam"]
  out$scale = scale
  out$criterion = criterion
  if (optModel) {
    kernel_list = list(type = kernel, par = opt_param["kparam"])
    anova_K = make_anovaKernel(x, x, kernel = kernel_list)
    K = combine_kernel(anova_K, theta)
    opt_model = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = opt_param["lambda"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "sramsvm"
  return(out)
}

theta_step.sramsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$opt_param["kparam"]
  n_class = object$n_class
  x = object$x
  y = object$y
  gamma = object$gamma
  theta = object$theta

  valid_x = object$valid_x
  valid_y = object$valid_y

  kernel_list = list(type = kernel, par = kparam)
  anova_K = make_anovaKernel(x, x, kernel = kernel_list)
  valid_anova_K = make_anovaKernel(valid_x, x, kernel_list)

  if (is.null(object$opt_model)) {
    K = combine_kernel(anova_K, theta)
    init_model = ramsvm_compact(K = K, y = y, gamma = gamma, lambda = lambda, ...)
  } else {
    init_model = object$opt_model
  }

  fold_err = mclapply(1:length(lambda_theta_seq),
                      function(j) {
                        error = try({
                          theta = find_theta.sramsvm(y = y, anova_kernel = anova_K, gamma = gamma, cmat = init_model$beta, c0vec = init_model$beta0,
                                                    n_class = n_class, lambda = lambda, lambda_theta = lambda_theta_seq[j], ...)
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

  out$opt_lambda_theta = opt_lambda_theta
  out$opt_ind = opt_ind
  out$opt_theta = opt_theta
  out$theta_seq = theta_seq
  out$opt_valid_err = opt_valid_err
  out$valid_err = valid_err

  if (isCombined) {
    optK = combine_kernel(anova_K, opt_theta)
    opt_model = ramsvm_compact(K = optK, y = y, gamma = gamma, lambda = lambda, ...)

  } else {
    opt_model = init_model
  }
  out$opt_model = opt_model
  class(out) = "sramsvm"
  return(out)
}


find_theta.sramsvm = function(y, anova_kernel, gamma, cmat, c0vec, n_class, lambda, lambda_theta, eig_tol_D = 0, epsilon_D = 1e-8)
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

  n = length(y)
  y_index = cbind(1:n, y)

  c0vec = as.matrix(c0vec)
  # convert y into ramsvm class code
  trans_Y = Y_matrix_gen(n_class, nobs = n, y = y)

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
    if(d == 1) {
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

































