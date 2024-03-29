mlapsvm_compact = function(K, L, y, lambda, lambda_I, epsilon = 1e-6,
                           eig_tol_D = 0, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-6, epsilon_I = 0)
{

  # The sample size, the number of classes and dimension of QP problem
  out = list()

  y_temp = factor(y)
  levs = levels(y_temp)
  attr(levs, "type") = class(y)
  y_int = as.integer(y_temp)

  n_class = length(levs)
  # if (is(y, "numeric")) {levs = as.numeric(classname)}

  n_l = length(y_int)
  n = nrow(K)
  n_u = n - n_l
  qp_dim = n_l * n_class

  # Convert y into msvm class code
  trans_Y = class_code(y_int, n_class)

  # Optimize alpha by solve.QP:
  # min (-d^Tb + 1/2 b^TDb)
  # subject to A^Tb <= b_0
  # Following steps (1) - (6)
  # (1) preliminary quantities
  J = cbind(diag(n_l), matrix(0, n_l, n - n_l))
  Jk = matrix(1, nrow = n_class, ncol = n_class)
  Ik = diag(1, n_class)

  # Vectorize y matrix
  y_vec = as.vector(trans_Y)

  # Index for non-trivial alphas
  nonzeroIndex = (y_vec != 1)

  # inv_LK = solve(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K))
  # LK = fixit(diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K), eig_tol_I)
  # inv_LK = chol2inv(chol(LK))
  # K = fixit(K, eig_tol_D)
  LK = diag(n_l * lambda, n) + n_l * lambda_I / n^2 * (L %*% K)
  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))
  # inv_LK = solve(LK + diag(max_LK * epsilon_I, n))
  # inv_LK = solve(LK + diag(max_LK * epsilon_I, n), t(J))
  # inv_LK = inverse(LK, epsilon = eig_tol_I)

  # inv_LK = solve(LK / max_LK + diag(epsilon_I, n), t(J) / max_LK)
  # inv_LK = solve(LK / max_LK + diag(epsilon_I, n), tol = eig_tol_I / 100) / max_LK
  inv_LK = solve(LK + diag(max(abs(LK)) * epsilon_I, n), tol = eig_tol_I)
  # inv_LK = chol2inv(chol(LK + diag(max_LK * epsilon_I, n)))

  # Q = J %*% K %*% inv_LK
  Q = J %*% K %*% inv_LK %*% t(J)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # Q = fixit(Q, epsilon = eig_tol_D)
  # Q = Q[1:n_l, 1:n_l]

  # (2) Compute D <- H
  D = (Ik - Jk / n_class) %x% Q

  # Subset the columns and rows for non-trivial alpha's

  Reduced_D = D[nonzeroIndex, nonzeroIndex]
  Reduced_D = fixit(Reduced_D, epsilon = eig_tol_D)
  max_D = max(abs(Reduced_D))
  # Reduced_D = Reduced_D / max_D
  diag(Reduced_D) = diag(Reduced_D) + max_D * epsilon_D

  # (3) Compute d <- g
  # g = -y_vec
  g = -y_vec

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

  # R consists of right hand sides of equality and inequality constraints
  r = c(r1, r2)

  # (6) Find solution by solve.QP.compact
  nonzero = find_nonzero(R)
  Amat = nonzero$Amat_compact
  Aind = nonzero$Aind
  dual = solve.QP.compact(Reduced_D, Reduced_g, Amat, Aind, r, meq = nrow(Reduced_R1))

  # Place the dual solution into the non-trivial alpha positions
  alpha = rep(0, qp_dim)
  alpha[nonzeroIndex] = dual$solution

  # Make alpha zero if they are too small
  alpha[alpha < 0] = 0
  # alpha[alpha > 1] = 1

  # Reshape alpha into a n by n_class matrix
  alpha = matrix(alpha, nrow = n_l)

  # Compute cmat = matrix of estimated coefficients

  cmat = -inv_LK %*% t(J) %*% (alpha - matrix(rep(rowMeans(alpha), n_class), ncol = n_class))
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

predict.mlapsvm_compact = function(object, newK = NULL)
{
  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels
  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}


mlapsvm = function(x = NULL, y, ux = NULL, lambda, lambda_I, kernel, kparam, scale = FALSE, adjacency_k = 6, normalized = FALSE,
                   weight = NULL, weightType = "Binary", epsilon = 1e-6,
                   eig_tol_D = 0, eig_tol_I = .Machine$double.eps, epsilon_D = 1e-8, epsilon_I = 0)
{
  out = list()
  n_l = NROW(x)
  n_u = NROW(ux)
  n = n_l + n_u
  rx = rbind(x, ux)
  p = ncol(x)
  center = rep(0, p)
  scaled = rep(1, p)

  if (scale) {
    rx = scale(rx)
    center = attr(rx, "scaled:center")
    scaled = attr(rx, "scaled:scale")
    x = (x - matrix(center, nrow = n_l, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_l, ncol = p, byrow = TRUE)
    ux = (ux - matrix(center, nrow = n_u, ncol = p, byrow = TRUE)) / matrix(scaled, nrow = n_u, ncol = p, byrow = TRUE)
  }

  K = kernelMatrix(rx, rx, kernel = kernel, kparam = kparam)

  # W = RSSL:::adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # d = rowSums(W)
  # L = diag(d) - W
  # if (normalized) {
  #   L = diag(1 / sqrt(d)) %*% L %*% diag(1 / sqrt(d))
  # }



  # W = adjacency_knn(rx, distance = "euclidean", k = adjacency_k)
  # graph = W
  graph = make_knn_graph_mat(rx, k = adjacency_k)
  L = make_L_mat(rx, kernel = kernel, kparam = kparam, graph = graph, weightType = weightType, normalized = normalized)


  solutions = mlapsvm_compact(K = K, L = L, y = y, lambda = lambda, lambda_I = lambda_I, epsilon = epsilon,
                              eig_tol_D = eig_tol_D, eig_tol_I = eig_tol_I, epsilon_D = epsilon_D, epsilon_I = epsilon_I)

  out$x = x
  out$ux = ux
  out$y = y
  out$n_class = solutions$n_class
  out$levels = solutions$levels
  out$weight = weight
  out$lambda = lambda
  out$lambda_I = lambda_I
  out$kparam = kparam
  out$cmat = solutions$cmat
  out$c0vec = solutions$c0vec
  out$alpha = solutions$alpha
  out$epsilon = epsilon
  out$eig_tol_D = eig_tol_D
  out$eig_tol_I = eig_tol_I
  out$epsilon_D = epsilon_D
  out$epsilon_I = epsilon_I
  out$kernel = kernel
  out$scale = scale
  out$center = center
  out$scaled = scaled
  out$fit_class = solutions$fit_class
  class(out) = "mlapsvm"
  return(out)
}


predict.mlapsvm = function(object, newx = NULL, newK = NULL)
{

  if (object$scale) {
    newx = (newx - matrix(object$center, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)) / matrix(object$scaled, nrow = nrow(newx), ncol = ncol(newx), byrow = TRUE)
  }

  if (is.null(newK)) {
    newK = kernelMatrix(newx, rbind(object$x, object$ux), kernel = object$kernel, kparam = object$kparam)
    # newK = kernelMatrix(rbfdot(sigma = object$kparam), newx, object$x)
  }

  cmat = object$cmat
  c0vec = object$c0vec
  levs = object$levels
  pred_y = (matrix(rep(c0vec, nrow(newK)), ncol = object$n_class, byrow = T) + (newK %*% cmat))
  pred_class = levs[apply(pred_y, 1, which.max)]

  if (attr(levs, "type") == "factor") {pred_class = factor(pred_class, levels = levs)}
  if (attr(levs, "type") == "numeric") {pred_class = as.numeric(pred_class)}
  if (attr(levs, "type") == "integer") {pred_class = as.integer(pred_class)}

  return(list(class = pred_class, pred_value = pred_y))
}

cv.mlapsvm = function(x, y, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5,
                      lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                      kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = 1, normalized = FALSE,
                      scale = FALSE, criterion = c("0-1", "loss"), optModel = FALSE, nCores = 1, ...)
{
  out = list()
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

  lambda_seq = as.numeric(lambda_seq)
  lambda_I_seq = as.numeric(lambda_I_seq)
  kparam = as.numeric(kparam)

  lambda_seq = sort(lambda_seq, decreasing = FALSE)
  lambda_I_seq = sort(lambda_I_seq, decreasing = TRUE)
  kparam = sort(kparam, decreasing = TRUE)

  # Combination of hyper-parameters
  # params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq[order(lambda_I_seq, decreasing = TRUE)], kparam = kparam)
  params = expand.grid(lambda = lambda_seq, lambda_I = lambda_I_seq, kparam = kparam)

  if (!is.null(valid_x) & !is.null(valid_y)) {
    model_list = vector("list", 1)
    fold_list = NULL

    #  Parallel computation on the combination of hyper-parameters
    fold_err = mclapply(1:nrow(params),
                        function(j) {
                          error = try({
                            msvm_fit = mlapsvm(x = x, y = y, ux = ux, lambda = params$lambda[j], lambda_I = params$lambda_I[j],
                                               kernel = kernel, kparam = params$kparam[j], scale = scale, normalized = normalized, ...)
                          })

                          if (!inherits(error, "try-error")) {
                            pred_val = predict.mlapsvm(msvm_fit, newx = valid_x)$class

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
    valid_err = round(sapply(fold_err, "[[", "error"), 8)
    model_list[[1]] = lapply(fold_err, "[[", "fit_model")
    opt_ind = max(which(valid_err == min(valid_err)))
    opt_param = params[opt_ind, ]
    opt_valid_err = min(valid_err)
  }

  out$opt_param = c(lambda = opt_param$lambda, lambda_I = opt_param$lambda_I, kparam = opt_param$kparam)
  out$opt_valid_err = opt_valid_err
  out$opt_ind = opt_ind
  out$valid_err = valid_err
  out$x = x
  out$y = y
  out$valid_x = valid_x
  out$valid_y = valid_y
  out$kernel = kernel
  out$scale = scale
  if (optModel) {
    opt_model = mlapsvm(x = x, y = y, ux = ux, lambda = opt_param$lambda, lambda_I = opt_param$lambda_I,
                        kernel = kernel, kparam = opt_param$kparam, scale = scale, normalized = normalized, ...)
    out$opt_model = opt_model
  }
  out$call = call
  class(out) = "mlapsvm"
  return(out)
}


