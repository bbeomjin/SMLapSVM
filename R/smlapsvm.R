# smlapsvm = function(x = NULL, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
#                     lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
#                     lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
#                     adjacency_k = 6, normalized = FALSE, weightType = "Binary",
#                     kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
#                     scale = TRUE, criterion = c("0-1", "loss"), isCombined = TRUE, nCores = 1, ...)
# {
#   out = list()
#   cat("Fit c-step \n")
#   cstep_fit = cstep.smlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
#                     lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = NULL,
#                     adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
#                     kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)
#
#   cat("Fit theta-step \n")
#   thetastep_fit = thetastep.smlapsvm(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)
#
#   cat("Fit c-step \n")
#   opt_cstep_fit = cstep.smlapsvm(x = x, y = y, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
#                     lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq, theta = thetastep_fit$opt_theta,
#                     adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
#                     kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = TRUE, nCores = nCores, ...)
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
#   out$ux = ux
#   out$n_class = opt_cstep_fit$n_class
#   class(out) = "smlapsvm"
#   return(out)
# }
#
# predict.smlapsvm = function(object, newx = NULL, newK = NULL)
# {
#   model = object$opt_model
#   cmat = model$cmat
#   c0vec = model$c0vec
#   levs = model$levels
#
#   # if (object$scale) {
#   #   newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
#   # }
#
#   if (is.null(newK)) {
#     new_anova_K = make_anovaKernel(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
#     newK = combine_kernel(new_anova_K, theta = object$opt_theta)
#     # newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
#     # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
#   }
#
#   pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = model$n_class, byrow = T) + (newK %*% cmat))
#   pred_class = levs[apply(pred_y, 1, which.max)]
#
#   if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
#   if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
#   if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}
#
#   return(list(class = pred_class, pred_value = pred_y))
# }



cstep.smlapsvm = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)}, theta = NULL,
                 adjacency_k = 6, normalized = FALSE, weightType = "Binary",
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

  if (is.null(theta)) {
    theta = rep(1, p)
  }

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  lambda_I_seq = sort(lambda_I_seq, decreasing = TRUE)
  kparam = sort(kparam, decreasing = FALSE)

  # Combination of hyper-parameters
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

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

    valid_err_mat = matrix(NA, nrow = length(kparam), ncol = nrow(params))

    for (i in 1:length(kparam)) {
      par = kparam[i]

      anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = par)
      # K = combine_kernel(anova_kernel = anova_K, theta = theta)

      # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
      # graph = W
      graph = make_knn_graph_mat(rx, k = adjacency_k)
      L = make_L_mat(rx, kernel = kernel, kparam = par, graph = graph, weightType = weightType, normalized = normalized)
      # L = fixit(L, epsilon = 0)
      #     if (any(theta > 0)) {
      # 	   graph = make_knn_graph_mat(rx[, theta > 0, drop = FALSE], k = adjacency_k)
      # 	   L = make_L_mat(rx[, theta > 0, drop = FALSE], kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)
      # 	  } else {
      # 	   graph = make_knn_graph_mat(rx, k = adjacency_k)
      # 	   L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)
      # 	  }

      valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = par)
      valid_K = combine_kernel(anova_kernel = valid_anova_K, theta = theta)
      #  Parallel computation on the combination of hyper-parameters
      fold_err = mclapply(1:nrow(params),
                          function(j) {
                            error = try({
                              msvm_fit = smlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y,
                                                          lambda = params$lambda[j], lambda_I = params$lambda_I[j], ...)
                            })

                            if (!inherits(error, "try-error")) {
                              pred_val = predict.mlapsvm_compact(msvm_fit, newK = valid_K)$class
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
      # model_list[[1]] = lapply(fold_err, "[[", "fit_model")
      valid_err_mat[i, ] = valid_err
    }
    opt_ind = which(valid_err_mat == min(valid_err_mat), arr.ind = TRUE)
    opt_ind = opt_ind[order(opt_ind[, 1], opt_ind[, 2], decreasing = c(FALSE, TRUE))[1], ]
    opt_param = c(lambda = params[opt_ind[2], 1], lambda_I = params[opt_ind[2], 2], kparam = kparam[opt_ind[1]])
    opt_valid_err = min(valid_err_mat)
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
  out$valid_x = valid_x
  out$valid_y = valid_y
  # out$adjacency_k = adjacency_k
  # out$normalized = normalized
  # out$weightType = weightType
  # out$anova_K = anova_K
  # out$K = K
  # out$valid_anova_K = valid_anova_K
  # out$valid_K = valid_K
  out$kernel = kernel
  out$kparam = opt_param["kparam"]
  out$scale = scale
  out$criterion = criterion
  if (optModel) {
    anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = opt_param["kparam"])
    opt_model = smlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, lambda = opt_param["lambda"], lambda_I = opt_param["lambda_I"], ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "smlapsvm"
  return(out)
}


thetastep.smlapsvm = function(object, lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                              isCombined = TRUE, optModel = FALSE, nCores = 1, ...)
{
  call = match.call()
  out = list()
  lambda_theta_seq = sort(as.numeric(lambda_theta_seq), decreasing = FALSE)
  lambda = object$opt_param["lambda"]
  lambda_I = object$opt_param["lambda_I"]
  criterion = object$criterion
  kernel = object$kernel
  kparam = object$opt_param["kparam"]
  n_class = object$n_class
  x = object$x
  y = object$y
  theta = object$theta
  ux = object$ux
  rx = rbind(x, ux)
  # adjacency_k = object$adjacency_k
  # normalized = object$normalized
  # weightType = object$weightType
  valid_x = object$valid_x
  valid_y = object$valid_y
  L = object$L
  # anova_K = object$anova_K
  # K = object$K

  anova_K = make_anovaKernel(rx, rx, kernel = kernel, kparam = kparam)
  valid_anova_K = make_anovaKernel(valid_x, rx, kernel = kernel, kparam = kparam)

  if (is.null(object$opt_model)) {
    init_model = smlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, lambda = lambda, lambda_I = lambda_I, ...)
  } else {
    init_model = object$opt_model
  }

  fold_err = mclapply(1:length(lambda_theta_seq),
                      function(j) {
                        error = try({
                          theta = find_theta.smlapsvm(y = y, anova_kernel = anova_K, L = L, cmat = init_model$cmat, c0vec = init_model$c0vec,
                                                      lambda = lambda, lambda_I = lambda_I, lambda_theta = lambda_theta_seq[j], ...)
                          if (isCombined) {
                            init_model = smlapsvm_compact(anova_K = anova_K, L = L, theta = theta, y = y, lambda = lambda, lambda_I = lambda_I, ...)
                          }
                        })

                        if (!inherits(error, "try-error")) {
                          valid_subK = combine_kernel(valid_anova_K, theta)
                          pred_val = predict.mlapsvm_compact(init_model, newK = valid_subK)$class

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

  if (optModel) {
    # subK = combine_kernel(anova_K, opt_theta)
    opt_model = smlapsvm_compact(anova_K = anova_K, L = L, theta = opt_theta, y = y, lambda = lambda, lambda_I = lambda_I, ...)
    out$opt_model = opt_model
  }

  class(out) = "smlapsvm"
  return(out)
}

find_theta.smlapsvm = function(y, anova_kernel, L, cmat, c0vec, lambda, lambda_I, lambda_theta = 1,
                               eig_tol_D = 0, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-8, epsilon_I = 1e-11)
{
  if (lambda_theta <= 0) {
    theta = rep(1, anova_kernel$numK)
    return(theta)
  }

  anova_kernel_orig = anova_kernel
  anova_kernel$K = lapply(anova_kernel$K, function(x) {
    diag(x) = diag(x) + max(abs(x)) * epsilon_I
    return(x)
  })

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)
  n_class = length(levs)

  n = NROW(cmat)
  n_l = length(y_int)
  n_u = n - n_l

  Y = class_code(y_int, k = n_class)

  Dmat = numeric(anova_kernel$numK)
  dvec = numeric(anova_kernel$numK)
  A_mat = NULL

  for(j in 1:anova_kernel$numK) {
    temp_D = 0
    temp_d = 0
    temp_A = NULL
    for (q in 1:ncol(cmat)) {
      cvec = cmat[, q]
      KLK_temp = anova_kernel_orig$K[[j]] %*% L %*% anova_kernel_orig$K[[j]]
      diag(KLK_temp) = diag(KLK_temp) + max(abs(KLK_temp)) * epsilon_I
      temp_D = temp_D + n_l * lambda_I / (2 * n^2) * t(cvec) %*% KLK_temp %*% cvec
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
  diag(Dmat) = diag(Dmat) + max_D * epsilon_D
  # Dmat = fixit(Dmat, epsilon = eig_tol_D, is_diag = TRUE)

  # Dmat = Dmat / max_D


  dvec_temp = matrix(1, nrow = n_l, ncol = n_class)
  dvec_temp[cbind(1:n_l, y_int)] = 0
  # dvec_temp = as.vector(Y)
  # dvec_temp[dvec_temp == 1] = 0
  # dvec_temp[dvec_temp < 0] = 1
  dvec = c(dvec, as.vector(dvec_temp))
  dvec = dvec
  # dvec = dvec / max_D

  # solve QP
  # diag(Dmat) = diag(Dmat) + epsilon_D

  A_mat = cbind(-A_mat, diag(1, n_l * n_class))
  A_mat = rbind(A_mat, diag(1, ncol(A_mat)))
  A_theta = cbind(diag(-1, anova_kernel$numK), matrix(0, anova_kernel$numK, (ncol(A_mat) - anova_kernel$numK)))
  A_mat = rbind(A_mat, A_theta)
  #    print(ncol(A.mat))
  bvec = c(rep(c0vec, each = n_l) - as.vector(Y), rep(0, anova_kernel$numK + n_l * n_class), rep(-1, anova_kernel$numK))
  #    print(A.mat)
  #    print(bvec)
  theta_sol = solve.QP(Dmat, -dvec, t(A_mat), bvec, meq = 0, factorized = FALSE)$solution
  theta = theta_sol[1:anova_kernel$numK]
  theta[theta < 1e-6] = 0
  theta = round(theta, 6)
  # theta_sol[theta_sol < 1e-6] = 0
  #    print(beta)
  return(theta)
}


smlapsvm_compact = function(anova_K, L, theta, y, lambda, lambda_I, epsilon = 1e-6,
                            eig_tol_D = 0, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-8, epsilon_I = 1e-11)
{

  # The sample size, the number of classes and dimension of QP problem
  out = list()

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)

  anova_K_orig = anova_K
  anova_K$K = lapply(anova_K$K, function(x) {
    diag(x) = diag(x) + max(abs(x)) * epsilon_I
    return(x)
  })

  K = combine_kernel(anova_K, theta = theta)

  if (sum(K) == 0) {
    diag(K) = 1
  }

  n = nrow(K)
  n_l = length(y_int)
  n_u = n - n_l
  qp_dim = n_l * n_class

  J = cbind(diag(n_l), matrix(0, n_l, n - n_l))

  KLK = 0
  for (i in 1:anova_K$numK) {
    KLK_temp = anova_K_orig$K[[i]] %*% L %*% anova_K_orig$K[[i]]
    diag(KLK_temp) = diag(KLK_temp) + max(abs(KLK_temp)) * epsilon_I
    KLK = KLK + theta[i]^2 * KLK_temp
  }

  lambda_K = n_l * lambda * K
  lambda_KLK = n_l * lambda_I / n^2 * KLK

  K_KLK = lambda_K + lambda_KLK
  # K_KLK = (K_KLK + t(K_KLK)) / 2

  inv_K_KLK = solve(K_KLK, tol = eig_tol_I)
  inv_K_KLK = (inv_K_KLK + t(inv_K_KLK)) / 2
  inv_K_KLK = inv_K_KLK %*% K %*% t(J)

  Q = J %*% K %*% inv_K_KLK
  # Q = fixit(Q, epsilon = eig_tol_D)
  # diag(Q) = diag(Q) + epsilon_D


  # Q = J %*% Q %*% t(J)
  # diag(Q) = diag(Q) + epsilon_D
  # Convert y into msvm class code
  trans_Y = class_code(y_int, n_class)

  # Optimize alpha by solve.QP:
  # min (-d^Tb + 1/2 b^TDb)
  # subject to A^Tb <= b_0
  # Following steps (1) - (6)
  # (1) preliminary quantities
  Jk = matrix(1, nrow = n_class, ncol = n_class)
  Ik = diag(1, n_class)

  # Vectorize y matrix
  y_vec = as.vector(trans_Y)

  # Index for non-trivial alphas
  nonzeroIndex = (y_vec != 1)

  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  # Q = K %*% inv_LK


  # Q = K %*% inv_KL

  # Q = Q[1:n_l, 1:n_l]

  # (2) Compute D <- H
  D = (Ik - Jk / n_class) %x% Q

  # Subset the columns and rows for non-trivial alpha's
  Reduced_D = D[nonzeroIndex, nonzeroIndex]
  Reduced_D = fixit(Reduced_D, epsilon = eig_tol_D)
  max_D = max(abs(Reduced_D))
  # Reduced_D = Reduced_D / max_D
  # diag(Reduced_D) = diag(Reduced_D) + epsilon_D
  diag(Reduced_D) = diag(Reduced_D) + max_D * epsilon_D

  # Reduced_D = nearPD(Reduced_D, eig.tol = rel_eig_tol)$mat
  # diag(Reduced_D) = diag(Reduced_D) + epsilon_D

  # (3) Compute d <- g
  g = -y_vec
  # g = -y_vec / max_D

  # Subset the components with non-trivial alpha's
  Reduced_g = g[nonzeroIndex]
  n_nonzeroIndex = length(Reduced_g)

  # (4) Compute A <- R
  # Equality constraint matrix
  R1 = ((Ik - Jk / n_class) %x% matrix(rep(1, n_l), nrow = 1))

  # Eliminate one redundant equality constraint
  R1 = matrix(R1[1:(n_class - 1), ], nrow = n_class - 1, ncol = ncol(R1))

  # Choose components with non-trivial alpha's
  Reduced_R1 = matrix(R1[, nonzeroIndex], nrow = nrow(R1), ncol = n_nonzeroIndex)

  # Inequality constraint matrix
  R2 = diag(rep(1, n_l * (n_class - 1)))
  R2 = rbind(R2, -R2)

  # R consists of equality and inequality constraints
  R = t(rbind(Reduced_R1, R2))

  # (5) compute (b_0, b) = r
  # Right hand side of equality constraints
  r1 = rep(0, nrow(Reduced_R1))

  # Right hand side of inequality constraints
  r2 = c(rep(0, nrow(R2) / 2), rep(-1, nrow(R2) / 2))
  # r2 = c(rep(0, nrow(R2) / 2), rep(-1 / n_l, nrow(R2) / 2))

  # R consists of right hand sides of equality and inequality constraints
  r = c(r1, r2)

  # (6) Find solution by solve.QP.compact
  nonzero = find_nonzero(R)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind

  dual = solve.QP.compact(Reduced_D, Reduced_g, Amat, Aind, r, meq = nrow(Reduced_R1))
  # dual = solve.QP(Reduced_D, Reduced_g, Amat, r, meq = nrow(Reduced_R1), factorized = TRUE)
  # dual_temp = solve.QP(Reduced_D, Reduced_g, R, r, meq = nrow(Reduced_R1))

  # Place the dual solution into the non-trivial alpha positions
  alpha = rep(0, qp_dim)
  alpha[nonzeroIndex] = dual$solution

  # Make alpha zero if they are too small
  alpha[alpha < 0] = 0
  # alpha[alpha > 1] = 1

  # Reshape alpha into a n by n_class matrix
  alpha = matrix(alpha, nrow = n_l)

  # Compute cmat = matrix of estimated coefficients

  # cmat = -inv_KLK %*% K %*% t(J) %*% (alpha - matrix(rep(rowMeans(alpha), n_class), ncol = n_class))
  cmat = -inv_K_KLK %*% (alpha - matrix(rep(rowMeans(alpha), n_class), ncol = n_class))
  # J = cbind(diag(1, n_l), matrix(0, n_l, n - n_l))

  # Find b vector
  Kcmat = J %*% K %*% cmat
  flag = T
  # while (flag) {
  #   logic = ((alpha > epsilon) & (alpha < (1 - epsilon)))
  #   # logic = ((alpha > epsilon) & (alpha < (1 - epsilon)))
  #   c0vec = numeric(n_class)
  #   if (all(colSums(logic) > 0)) {
  #     # Using alphas between 0 and 1, we get c0vec by KKT conditions
  #     for (i in 1:n_class) {
  #       c0vec[i] = mean((trans_Y[, i] - Kcmat[, i])[logic[, i]])
  #     }
  #     if (abs(sum(c0vec)) < 0.001) {
  #       flag = F
  #     } else {
  #       epsilon = min(epsilon * 2, 0.5)
  #     }
  #   } else {
      flag = F
      # Otherwise, LP starts to find b vector
      # reformulate LP w/o equality constraint and redudancy
      # objective function with (b_j)_+,-(b_j)_, j=1,...,(k-1) and \xi_ij

      a = c(rep(0, 2 * (n_class - 1)), rep(1, n_l * (n_class - 1)))
      # inequality conditions
      B1 = diag(-1, n_class - 1) %x% rep(1, n_l)
      B2 = matrix(1, n_l, n_class - 1)
      A = cbind(B1, -B1)
      A = rbind(A, cbind(B2, -B2))
      A = A[nonzeroIndex, ] # reduced.A
      A = cbind(A, diag(1, n_l * (n_class - 1)))
      b = matrix(Kcmat - trans_Y, ncol = 1)
      b = b[nonzeroIndex] # reduced.b
      # constraint directions
      const.dir = matrix(rep(">=", nrow(A)))

      bpos = lp("min", objective.in = a, const.mat = A, const.dir = const.dir,
                const.rhs = b)$solution[1:(2 * (n_class - 1))]
      c0vec = cbind(diag(1, n_class - 1), diag(-1, n_class - 1)) %*% matrix(bpos, ncol = 1)
      c0vec = c(c0vec, -sum(c0vec))
    # }
  # }

  # Compute the fitted values
  fit = (matrix(rep(c0vec, n_l), ncol = n_class, byrow = T) + Kcmat)
  fit_class = levs[apply(fit, 1, which.max)]
  if (attr(levs, "type") == "factor") {fit_class = factor(fit_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {fit_class = as.numeric(fit_class)}
  if (attr(levs, "type") == "integer") {fit_class = as.integer(fit_class)}
  # table(y, fit_class)
  # Return the output
  out$alpha = alpha
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








