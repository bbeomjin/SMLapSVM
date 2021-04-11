srmlapsvm = function(x = NULL, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5,
                    lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                    adjacency_k = 6, normalized = TRUE, weightType = "Binary",
                    kernel = c("linear", "radial", "poly", "spline", "anova_radial"), kparam = c(1),
                    scale = FALSE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, ...)
{
  out = list()
  cat("Fit c-step \n")
  cstep_fit = cstep.rmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
                    adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                    kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
  
  cat("Fit theta-step \n")
  theta_step_fit = theta_step.srmlapsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)
  
  cat("Fit c-step \n")
  opt_cstep_fit = cstep.rmlapsvm(x = x, y = y, ux = ux, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                        lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = theta_step_fit$opt_theta,
                        adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
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
    newK = combine_kernel(new_anova_K, theta = object$theta)
    # newK = kernelMat(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }
  
  cmat = model$cmat
  c0vec = model$c0vec
  
  pred_y = (matrix(c0vec, nrow = nrow(newK), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = apply(pred_y, 1, which.max)
  return(list(class = pred_class, pred_value = pred_y))
}


cstep.rmlapsvm = function(x, y, ux = NULL, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5, 
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
    K = combine_kernel(anova_kernel = anova_K, theta = theta)
    
    W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
    # graph = make_knn_graph_mat(rx, k = adjacency_k)
    graph = W
    L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)
    
    valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel_list)
    valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)
    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          msvm_fit = rmlapsvm_compact(K = K, L = L, y = y, gamma = gamma, 
                                                      lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                          
                          pred_val = predict.rmlapsvm_compact(msvm_fit, newK = valid_K)$class
                          
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
  out$gamma = gamma
  out$L = L
  out$n_class = n_class
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$anova_K = anova_K
  out$K = K
  out$valid_anova_K = valid_anova_K
  out$valid_K = valid_K
  out$kernel = kernel
  out$scale = scale
  out$criterion = criterion
  if (optModel) {
    opt_model = rmlapsvm_compact(K = K, L = L, y = y, gamma = gamma, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, ...)
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
  # ux = object$ux
  # rx = rbind(x, ux)
  valid_y = object$valid_y
  
  anova_K = object$anova_K
  K = object$K
  L = object$L
  valid_anova_K = object$valid_anova_K
  if (is.null(object$opt_model)) {
    init_model = rmlapsvm_compact(K = K, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
  } else {
    init_model = object$opt_model
  }
  
  fold_err = mclapply(1:length(lambda_theta_seq),
                      function(j) {
                        theta = find_theta.srmlapsvm(y = y, gamma = gamma, anova_kernel = anova_K, L = L, 
                                                     cmat = init_model$cmat, c0vec = init_model$c0vec, n_class = n_class, 
                                                     lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j])
                        
                        if (isCombined) {
                          subK = combine_kernel(anova_K, theta)
                          init_model = rmlapsvm_compact(K = subK, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
                        }
                        
                        valid_subK = combine_kernel(valid_anova_K, theta)
                        pred_val = predict.rmlapsvm_compact(init_model, newK = valid_subK)$class
                        
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
    subK = combine_kernel(anova_K, opt_theta)
    opt_model = rmlapsvm_compact(subK, L = L, y = y, gamma = gamma, lambda = lambda, lambda_I = lambda_I, ...)
    
  } else {
    opt_model = init_model
  }
  out$opt_model = opt_model
  return(out)
}


find_theta.srmlapsvm = function(y, gamma, anova_kernel, L, cmat, c0vec, n_class, lambda, lambda_I, lambda_theta = 1, epsilon = 1e-6)
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
  Dmat = diag(Dmat)
  
  # dvec_temp = matrix(1, nrow = n_l, ncol = n_class)
  # dvec_temp[cbind(1:n_l, y)] = 0
  
  dvec_temp = matrix(1 - gamma, nrow = n_l, ncol = n_class)
  dvec_temp[y_index] = gamma
  
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




















































