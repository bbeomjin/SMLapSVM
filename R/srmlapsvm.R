srmlapsvm = function(x = NULL, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                    lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                    adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                    kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                    scale = FALSE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, verbose = 0, ...)
{
  out = list()

  fun_ind = ifelse(gamma == 0.0, 1, 2)

  cstep_fun = switch(fun_ind,
                     cstep.smlapsvm,
                     cstep.srmlapsvm)

  theta_step_fun = switch(fun_ind,
                          theta_step.smlapsvm,
                          theta_step.srmlapsvm)


  cat("Fit c-step \n")
  # cstep_fit = cstep.srmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
  #                   lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
  #                   adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
  #                   kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  cstep_args = list(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
                    adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                    kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  if (fun_ind == 2) {cstep_args$gamma = gamma}
  cstep_fit = do.call(cstep_fun, cstep_args)

  cat("Fit theta-step \n")
  theta_step_fit = theta_step_fun(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)


  cat("Fit c-step \n")
  opt_cstep_args = list(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                        lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = theta_step_fit$opt_theta,
                        adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                        kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
  # opt_cstep_fit = cstep.srmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
  #                       lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = theta_step_fit$opt_theta,
  #                       adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
  #                       kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
  if (fun_ind == 2) {opt_cstep_args$gamma = gamma}
  opt_cstep_fit = do.call(cstep_fun, opt_cstep_args)

  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
	cat("CV-error(theta-step):", theta_step_fit$opt_valid_err, "\n")
	cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
  }

  out$opt_param = opt_cstep_fit$opt_param
  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param$kparam
  out$opt_theta = theta_step_fit$opt_theta
  out$theta = theta_step_fit$theta
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
  kernel = object$kernel
  kparam = object$kparam
  cmat = model$cmat
  c0vec = model$c0vec

  # if (object$scale) {
  #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  # }

  if (is.null(newK)) {
    new_anova_K = make_anovaKernel(newx, rbind(object$x, object$ux), kernel = list(type = object$kernel, par = object$kparam))
    newK = combine_kernel(new_anova_K, theta = object$opt_theta)
    # newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = model$cmat
  c0vec = model$c0vec

  pred_y = (matrix(c0vec, nrow = nrow(newK), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = apply(pred_y, 1, which.max)
  return(list(class = pred_class, pred_value = pred_y))
}


cstep.srmlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)}, theta = NULL,
                 adjacency_k = 6, normalized = FALSE, weightType = "Binary",
                 kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                 scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  out = list()
  p = ncol(x)

  if (!is.numeric(lambda_seq)) {
    lambda_seq = as.numeric(lambda_seq)
  }

  if (!is.numeric(lambda_I_seq)) {
    lambda_I_seq = as.numeric(lambda_I_seq)
  }

  if (!is.numeric(kparam)) {
    kparam = as.numeric(kparam)
  }

  if (is.null(theta)) {
    theta = rep(1, p)
  }


  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq, kparam = kparam)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    n_l = NROW(x)
    n_u = NROW(ux)
    n = n_l + n_u
    rx = rbind(x, ux)


    # The number of classes
    n_class = length(unique(y))

    center = rep(0, p)
    scaled = rep(1, p)

    if (scale) {
      rx = scale(rx)
      center = attr(rx, "scaled:center")
      scaled = attr(rx, "scaled:scale")
      x = (x - matrix(center, nrow = n_l, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_l, ncol = p, byrow = TRUE)
      ux = (ux - matrix(center, nrow = n_u, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_u, ncol = p, byrow = TRUE)
    }

    kernel_list = list(type = kernel, par = kparam)
    anova_K = make_anovaKernel(rx, rx, kernel = kernel_list)
    # K = combine_kernel(anova_kernel = anova_K, theta = theta)

    # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
    # graph = W

	  graph = make_knn_graph_mat(rx, k = adjacency_k)
    L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)

    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel_list)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)
    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            msvm_fit = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                         lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.rmlapsvm_compact(msvm_fit, newK = valid_K)$class
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
  out$opt_param = opt_param
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$ux = ux
  out$y = y
  out$theta = theta
  out$gamma = gamma
  out$L = L
  out$n_class = n_class
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$anova_K = anova_K
  # out$K = K
  out$valid_anova_K = valid_anova_K
  out$valid_K = valid_K
  out$kernel = kernel
  out$scale = scale
  out$criterion = criterion
  if (optModel) {
    opt_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                  lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, ...)
    out$opt_model = opt_model
  }
  out$call = call
  return(out)
}


theta_step.srmlapsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda = object$opt_param$lambda
  lambda_I = object$opt_param$lambda_I
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$opt_param$kparam
  n_class = object$n_class
  gamma = object$gamma
  # x = object$x
  y = object$y
  theta = object$theta
  # ux = object$ux
  # rx = rbind(x, ux)
  valid_y = object$valid_y

  anova_K = object$anova_K
  K = object$K
  L = object$L
  valid_anova_K = object$valid_anova_K
  if (is.null(object$opt_model)) {
    init_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
  } else {
    init_model = object$opt_model
  }

  fold_err = mclapply(1:length(lambda_theta_seq),
                      function(j) {
                        error = try({
                          theta = find_theta.srmlapsvm(y = y, gamma = gamma, anova_kernel = anova_K, L = L,
                                                       cmat = init_model$cmat, c0vec = init_model$c0vec, n_class = n_class,
                                                       lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
                          if (isCombined) {
                            # subK = combine_kernel(anova_K, theta)
                            init_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, gamma = gamma,
                                                           lambda = lambda, lambda_I = lambda_I, ...)
                          }
                        })

                        if (!inherits(error, "try-error")) {
                          valid_subK = combine_kernel(valid_anova_K, theta)
                          pred_val = predict.rmlapsvm_compact(init_model, newK = valid_subK)$class

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
    # subK = combine_kernel(anova_K, opt_theta)
    opt_model = srmlapsvm_compact(anova_K = anova_K, L = L, theta = opt_theta, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)

  } else {
    opt_model = init_model
  }
  out$opt_model = opt_model
  return(out)
}


find_theta.srmlapsvm = function(y, gamma, anova_kernel, L, cmat, c0vec, n_class, lambda, lambda_I, lambda_theta = 1,
                                eig_tol_D = .Machine$double.eps, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-13)
{
  n = NROW(cmat)
  n_l = length(y)
  n_u = n - n_l

  # Y = class_code(y, k = n_class)
  y_index = cbind(1:n_l, y)

  Dmat = numeric(anova_kernel$numK)
  dvec = numeric(anova_kernel$numK)
  A_mat = NULL

  for(j in 1:anova_kernel$numK) {
    temp_D = 0
    temp_d = 0
    temp_A = NULL
    for (q in 1:ncol(cmat)) {
      cvec = cmat[, q]
      temp_D = temp_D + n_l * lambda_I / (n_l + n_u)^2 * t(cvec) %*% anova_kernel$K[[j]] %*% L %*% anova_kernel$K[[j]] %*% cvec
      temp_d = temp_d + n_l * lambda / 2 * t(cvec) %*% anova_kernel$K[[j]] %*% cvec + n_l * lambda_theta
      temp_A = rbind(temp_A, (anova_kernel$K[[j]][1:n_l, ] %*% cvec))
    }
    Dmat[j] = temp_D
    dvec[j] = temp_d
    A_mat = cbind(A_mat, temp_A)
  }

  Dmat = c(Dmat, c(rep(0, n_l * n_class)))
  # max_D = max(abs(Dmat))
  Dmat = diag(Dmat)
  Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)
  # Dmat = nearPD(Dmat)$mat
  # Dmat = Dmat / max_D

  # dvec_temp = matrix(1, nrow = n_l, ncol = n_class)
  # dvec_temp[cbind(1:n_l, y)] = 0

  dvec_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  dvec_temp[y_index] = gamma

  # dvec_temp = as.vector(Y)
  # dvec_temp[dvec_temp == 1] = 0
  # dvec_temp[dvec_temp < 0] = 1
  dvec = c(dvec, as.vector(dvec_temp))
  dvec = dvec
  # dvec = dvec / max_D

  # solve QP

  # diag(Dmat) = diag(Dmat) + epsilon_D

  m_index = matrix(1:(n_l * n_class), ncol = n_class)[cbind(1:n_l, y)]
  A_mat[m_index, ] = -A_mat[m_index, ]
  A_mat = cbind(-A_mat, diag(1, n_l * n_class))
  A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)

  bb = c0vec[y]
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(c0vec, nrow = n_l, ncol = n_class, byrow = TRUE)
  bb_j[y_index] = bb_yi
  bvec = c(as.vector(bb_j), rep(0, anova_kernel$numK + n_l * n_class), rep(-1, anova_kernel$numK))

  theta_sol = solve.QP(Dmat, -dvec, t(A_mat), bvec, meq = 0, factorized = FALSE)$solution
  theta = theta_sol[1:anova_kernel$numK]
  theta[theta < 1e-8] = 0

  return(theta)
}



srmlapsvm_compact = function(anova_K, L, theta, y, gamma = 0.5, lambda, lambda_I,
                             epsilon = 1e-6, eig_tol_D = .Machine$double.eps, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-13)
{
  out = list()
  # The labeled sample size, unlabeled sample size, the number of classes and dimension of QP problem
  n_class = length(unique(y))

  K = combine_kernel(anova_K, theta = theta)
  # K = (K + t(K)) / 2

  if (sum(K) == 0) {
    diag(K) = 1
  }

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

  m_mat = 0
  for (i in 1:anova_K$numK) {
    m_mat = m_mat + n_l * lambda_I / n^2 * theta[i]^2 * anova_K$K[[i]] %*% L %*% anova_K$K[[i]]
  }

  # K = fixit(K, eig_tol_D)
  # m_mat = fixit(m_mat, eig_tol_D)

  KLK = n_l * lambda * K + m_mat
  # KLK = (KLK + t(KLK)) / 2
  # KLK = fixit(KLK, epsilon = eig_tol_I)
  # KLK = nearPD(KLK, eig.tol = rel_eig_tol)$mat
  # diag(KLK) = diag(KLK) + epsilon_D
  # inv_KLK = solve(KLK)
  # inv_KLK = chol2inv(chol(KLK))
  inv_KLK = inverse(KLK, epsilon = eig_tol_I)

  # inv_KLK = solve(n_l * lambda * K + m_mat + diag(epsilon, n))

  Q = J %*% K %*% inv_KLK %*% K %*% t(J)
  Q = fixit(Q, epsilon = 0)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # diag(Q) = diag(Q) + epsilon_D

  # Compute Q = K x inv_LK
  D = matrix(0, qp_dim, qp_dim)
  Amat = matrix(0, (2 * qp_dim + n_class), qp_dim)

  for (k in 1:n_class) {
    D = D + t(Hmatj[[k]]) %*% Q %*% Hmatj[[k]]
    Amat[k, ] = rep(1, n_l) %*% Hmatj[[k]]
  }
  # D = fixit(D)
  # D = fixit(D, epsilon = eig_tol_D)
  max_D = max(abs(D))
  D = D / max_D
  diag(D) = diag(D) + epsilon_D

  # D = nearPD(D, eig.tol = rel_eig_tol)$mat
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

  # dvec = -g
  dvec = -g / max_D

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
  #       bvec[n_class + qp_dim + (j - 1) * n_l + i] = bvec[n_class + qp_dim + (j - 1) * n_l + i] - epsilon
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
  # alpha_mat[y_index][alpha_mat[y_index] > gamma] = gamma
  #
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

  cmat_temp = matrix(0, n_l, n_class)
  for (k in 1:n_class) {
    cmat_temp[, k] = Hmatj[[k]] %*% alpha
  }
  cmat = inv_KLK %*% K %*% t(J) %*% cmat_temp

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

  ######################### 수정필요 ########################################
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

  # print(dim(Alp))
  # print(length(blp))


  # Alp = matrix(0, nrow = qp_dim + 1, ncol = (qp_dim + 2 * n_class))
  # blp = rep(0, qp_dim + 1)
  #
  # for (j in 1:n_class) {
  #   Alp[1, (qp_dim + 2 * j - 1)] = 1
  #   Alp[1, (qp_dim + 2 * j)] = -1
  # }
  #
  # for(j in 1:n_class) {
  #   for(i in 1:n_l) {
  #     Alp[(1 + n_l * (j - 1) + i), n_l * (j - 1) + i] = 1
  #     if (y[i] == j) {
  #       Alp[(1 + n_l * (j - 1) + i), (qp_dim + 2 * (j - 1) + 1)] = 1
  #       Alp[(1 + n_l * (j - 1) + i), (qp_dim + 2 * (j - 1) + 2)] = -1
  #       blp[(1 + n_l * (j - 1) + i)] = (k - 1) - Kcmat[i, j]
  #     }
  #     if (y[i] != j) {
  #       Alp[(1 + n_l * (j - 1) + i), (qp_dim + 2 * (j - 1) + 1)] = -1
  #       Alp[(1 + n_l * (j - 1) + i), (qp_dim + 2 * (j - 1) + 2)] = 1
  #       blp[(1 + n_l * (j - 1) + i)] = 1 + Kcmat[i, j]
  #     }
  #   }
  # }
  # print(dim(Alp))
  # print(length(blp))

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

















































