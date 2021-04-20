sramlapsvm = function(x = NULL, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                    lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                    gamma = 0.5, adjacency_k = 6, normalized = FALSE, weightType = "Binary",
                    kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                    scale = TRUE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep.sramlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
                    gamma = gamma, adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                    kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  cat("Fit theta-step \n")
  theta_step_fit = theta_step.sramlapsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  cat("Fit c-step \n")
  opt_cstep_fit = cstep.sramlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                        lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = theta_step_fit$opt_theta,
                        gamma = gamma, adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                        kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)

  out$opt_valid_err = opt_cstep_fit$opt_valid_err
  out$opt_model = opt_cstep_fit$opt_model
  out$kernel = kernel
  out$kparam = opt_cstep_fit$opt_param$kparam
  out$theta = theta_step_fit$opt_theta
  out$x = x
  out$y = y
  out$ux = ux
  out$n_class = opt_cstep_fit$n_class
  class(out) = "sramlapsvm"
  return(out)
}

predict.sramlapsvm = function(object, newx = NULL, newK = NULL)
{
  model = object$opt_model
  kernel = object$kernel
  kparam = object$kparam
  cmat = model$beta
  c0vec = model$beta0

  # if (object$scale) {
  #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  # }

  if (is.null(newK)) {
    new_anova_K = make_anovaKernel(newx, rbind(object$x, object$ux), kernel = list(type = object$kernel, par = object$kparam))
    newK = combine_kernel(new_anova_K, theta = object$theta)
    # newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  temp_pred_y = predict_kernel(K_test = newK,
                               beta = cmat,
                               beta0 = c0vec,
                               k = object$n_class)
  inner_prod = temp_pred_y$inner_prod
  pred_class = temp_pred_y$class

  return(list(class = pred_class, pred_value = inner_prod))
}


cstep.sramlapsvm = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)}, gamma = 0.5,
                 theta = NULL, adjacency_k = 6, normalized = FALSE, weightType = "Binary",
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

    # fit = angle_lapsvm_core(K = K, L = L, y = y, lambda = params$lambda[j], lambda_I = params$lambda_I[j], maxiter = 1e+4)

    # predict.angle_lapsvm_core(fit, newK = K[1:30, ])[[2]]
    # table(predict.angle_lapsvm_core(fit, newK = K[1:30, ])[[1]], y)


    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          msvm_fit = sramlapsvm_core(anova_K = anova_K, L = L, theta = theta, y = y, lambda = params$lambda[j], lambda_I = params$lambda_I[j], gamma = gamma, ...)
                          # msvm_fit = angle_lapsvm_core(K = K, L = L, y = y, lambda = params$lambda[j], lambda_I = params$lambda_I[j], gamma = gamma)

                          pred_val = predict.ramlapsvm_core(msvm_fit, newK = valid_K)$class

                          if (criterion == "0-1") {
                            acc = sum(valid_y == pred_val) / length(valid_y)
                            err = 1 - acc
                          } else {
                            # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
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
  out$L = L
  out$theta = theta
  out$gamma = gamma
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
    opt_model = sramlapsvm_core(anova_K = anova_K, L = L, theta = theta, y = y, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, gamma = gamma, ...)
    # opt_model = angle_lapsvm_core(K = K, L = L, y = y, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, gamma = gamma)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "sramlapsvm"
  return(out)
}

theta_step.sramlapsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)}, isCombined = TRUE, nCores = 1, ...)
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
    init_model = sramlapsvm_core(anova_K = anova_K, L = L, theta = theta, y = y, lambda = lambda, lambda_I = lambda_I, ...)
  } else {
    init_model = object$opt_model
  }

  fold_err = mclapply(1:length(lambda_theta_seq),
                      function(j) {
                        theta = find_theta.sramlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = init_model$beta, c0vec = init_model$beta0,
                                                      gamma = gamma, n_class = n_class, lambda = lambda, lambda_I = lambda_I,
                                                      lambda_theta = lambda_theta_seq[j])

                        if (isCombined) {
                          # subK = combine_kernel(anova_K, theta)
                          init_model = sramlapsvm_core(anova_K = anova_K, L = L, theta = theta, y = y, lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)
                          # init_model = angle_lapsvm_core(K = subK, L = L, y = y, lambda = lambda, lambda_I = lambda_I, gamma = gamma)
                        }

                        valid_subK = combine_kernel(valid_anova_K, theta)
                        pred_val = predict.ramlapsvm_core(init_model, newK = valid_subK)$class

                        if (criterion == "0-1") {
                          acc = sum(valid_y == pred_val) / length(valid_y)
                          err = 1 - acc
                        } else {
                          # err = ramsvm_hinge(valid_y, pred_val$inner_prod, k = k, gamma = gamma)
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
    opt_model = sramlapsvm_core(anova_K = anova_K, L = L, theta = opt_theta, y = y, lambda = lambda, lambda_I = lambda_I, gamma = gamma, ...)

  } else {
    opt_model = init_model
  }
  out$opt_model = opt_model
  class(out) = "sramlapsvm"
  return(out)
}


find_theta.sramlapsvm = function(y, anova_kernel, L, cmat, c0vec, gamma, n_class, lambda, lambda_I, lambda_theta = 1, epsilon = 1e-6)
{

  if (anova_kernel$numK == 1)
  {
    cat("Only one kernel", "\n")
    return(c(1))
  }

  # Standard QP form :
  # min
  n = NROW(cmat)
  n_l = length(y)
  n_u = n - n_l
  n_class = length(unique(y))

  trans_Y = Y_matrix_gen(n_class, nobs = n_l, y = y)
  Y_code = XI_gen(n_class)

  cmat = as.matrix(cmat)
  c0vec = as.vector(c0vec)

  Dmat = numeric(anova_kernel$numK)
  dvec = numeric(anova_kernel$numK)
  A_mat = NULL

  for(j in 1:anova_kernel$numK) {
    temp_D = 0
    temp_d = 0
    # temp_A = NULL
    for (q in 1:ncol(cmat)) {
      cvec = cmat[, q]
      temp_D = temp_D + n_l * lambda_I / (n_l + n_u)^2 * t(cvec) %*% anova_kernel$K[[j]] %*% L %*% anova_kernel$K[[j]] %*% cvec
      temp_d = temp_d + n_l * lambda / 2 * t(cvec) %*% anova_kernel$K[[j]] %*% cvec + n_l * lambda_theta
    }
    temp_A = as.vector((anova_kernel$K[[j]][1:n_l, ] %*% cmat) %*% Y_code)

    Dmat[j] = temp_D
    dvec[j] = temp_d
    A_mat = cbind(A_mat, temp_A)
  }

  Dmat = c(Dmat, c(rep(0, n_l * n_class)))
  Dmat = diag(Dmat)

  dvec_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  dvec_temp[cbind(1:n_l, y)] = gamma
  # dvec_temp = as.vector(Y)
  # dvec_temp[dvec_temp == 1] = 0
  # dvec_temp[dvec_temp < 0] = 1
  dvec = c(dvec, as.vector(dvec_temp))

  # solve QP
  diag(Dmat) = diag(Dmat) + epsilon

  m_index = matrix(1:(n_l * n_class), ncol = n_class)[cbind(1:n_l, y)]
  A_mat[m_index, ] = -A_mat[m_index, ]
  A_mat = cbind(-A_mat, diag(1, n_l * n_class))
  A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)
  #    print(ncol(A.mat))

  bb = rowSums(sapply(1:(n_class - 1), function(x) trans_Y[, x] * c0vec[x]))
  bb_yi = (n_class - 1) - bb
  bb_j = 1 + matrix(rowSums(Y_code * matrix(c0vec, nrow = nrow(Y_code), ncol = ncol(Y_code), byrow = TRUE)), nrow = n_l, ncol = n_class, byrow = TRUE)
  bb_j[cbind(1:n_l, y)] = bb_yi
  bvec = c(as.vector(bb_j), rep(0, anova_kernel$numK + n_l * n_class), rep(-1, anova_kernel$numK))

  #    print(A.mat)
  #    print(bvec)
  theta_sol = solve.QP(Dmat, -dvec, t(A_mat), bvec, meq = 0, factorized = FALSE)$solution
  theta = theta_sol[1:anova_kernel$numK]
  theta[theta < 1e-8] = 0
  # theta_sol[theta_sol < 1e-6] = 0
  #    print(beta)
  return(theta)
}


sramlapsvm_core = function(anova_K, L, theta, y, gamma = 0.5, lambda, lambda_I, weight = NULL, epsilon = 1e-4 * length(y) * length(unique(y)), maxiter = 300, epsilon_D = 1e-6)
{

  out = list()
  n_class = length(unique(y))
  n_l = length(y)

  K = combine_kernel(anova_K, theta = theta)

  n = nrow(K)
  n_u = n - n_l

  m_mat = 0
  for (i in 1:anova_K$numK) {
    m_mat = m_mat + n_l * lambda_I / n^2 * theta[i]^2 * anova_K$K[[i]] %*% L %*% anova_K$K[[i]]
  }

  inv_KLK = solve(n_l * lambda * K + m_mat + diag(epsilon_D, n))

  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  Q = (n_l * lambda) * (K %*% inv_KLK %*% K)[1:n_l, 1:n_l] + 1

  warm = matrix(data = 0.0, nrow = n_l, ncol = n_class)
  if (is.null(weight)) {weight = numeric(n_l) + 1.0}

  #------------------------------------------------------------------#
  # Create k-vertex simplex.                                         #
  #------------------------------------------------------------------#
  my = t(XI_gen(k = n_class))

  yyi = Y_matrix_gen(k = n_class, nobs = n_l, y = y)

  alpha_ij = warm
  alpha_yi = numeric(n_l)


  erci = as.double(-rep(1, ncol(Q)) / n_l / lambda)

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
                     inv_LK = inv_KLK %*% K)
  # drop(crossprod(beta[[1]][, 1], X))

  # tt = beta_linear(x = x, y = y_train, k = k, my = my, warm = warm, lambda = templambda)

  # beta0 = matrix(find_theta2(y, K, gamma = gamma, cmat = beta$beta, lambda = lambda), nrow = 1)

  betaout = beta$beta
  beta0out = beta$beta0
  # beta0out[[count]] = beta0

  out$y = y
  out$n_class = n_class
  out$gamma = gamma
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$beta = betaout
  out$beta0 = beta0out
  out$epsilon = epsilon
  out$warm = warm
  class(out) = "ramlapsvm_core"
  return(out)
}




