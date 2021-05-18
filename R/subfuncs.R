main_kernel = function(x, y, kernel)
{
  x = as.matrix(x)
  y = as.matrix(y)
  if (kernel$type == "linear")
    K = (x %*% t(y))
  if (kernel$type == "poly")
    K = (1 + x %*% t(y))^kernel$par
  if (kernel$type == "radial" | kernel$type == "radial2")
  {
    # normx = drop((x^2) %*% rep(1.0, ncol(x)))
    normx = rowSums(x^2)
    # normy = drop((y^2) %*% rep(1.0, ncol(y)))
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kernel$par)
    # K_temp = kernlab::kernelMatrix(rbfdot(sigma = kernel$par), as.matrix(x), as.matrix(y))
  }
  # if (sym) {
  #   K = (K + t(K)) / 2
  # }
  return(K)
}

kernelMat = function(x, y, kernel = "radial", kparam = 1.0) {

  if (NCOL(x) == 0) {
    x = matrix(1, nrow = nrow(x), ncol = 1)
  }

  if (NCOL(y) == 0) {
    y = matrix(1, nrow = nrow(y), ncol = 1)
  }

  if (kernel == "poly") {
    obj = (x %*% t(y) + 1.0)^kparam
  } else if (kernel == "radial" | kernel == "radial2") {
    normx = rowSums(x^2)
    # normy = drop((y^2) %*% rep(1.0, ncol(y)))
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    obj = exp(-temp * kparam)
  } else if (kernel == "spline") {
    K = 0
    p = ncol(x)
    spline_kernel = function(x, u)
    {
      x = as.matrix(x)
      u = as.matrix(u)
      K1x = (x - 1/2)
      K1u = (u - 1/2)
      K2x = (K1x^2 - 1/12) / 2
      K2u = (K1u^2 - 1/12) / 2
      ax = x %x% matrix(1, 1, nrow(u))
      au = u %x% matrix(1, 1, nrow(x))
      b = abs(ax - t(au))
      K1 = K1x %x% t(K1u)
      K2 = K2x %x% t(K2u) - ((b - 1 / 2)^4 - (b - 1 / 2)^2 / 2 + 7 / 240) / 24
      return(list(K1 = K1, K2 = K2))
    }
    for(d in 1:p)
    {
      K_temp = spline_kernel(as.matrix(x[, d]), as.matrix(y[, d]))
      K = K + K_temp$K1 + K_temp$K2
    }
    obj = K
  } else if (kernel == "linear") {
    obj = tcrossprod(x, y)
  } else if (kernel == "anova_radial") {
    K = 0
    for (d in 1:NCOL(x))
    {
      A = as.matrix(x[, d])
      B = as.matrix(y[, d])
      K_temp = main_kernel(A, B, kernel = list(type = "radial", par = kparam))
      K = K + K_temp
    }
    obj = K
  } else {
    obj = NULL
  }
  # if (sym) {
  #   obj = (obj + t(obj)) / 2
  # }
  return(obj)
}

XI_gen = function(k) {

  tempA = -(1.0 + sqrt(k)) / ((k - 1)^(1.5))
  tempB = tempA + sqrt(k / (k - 1))

  XI = matrix(data = tempA, nrow = k - 1, ncol = k)

  XI[, 1] = 1.0 / sqrt(k - 1)

  for (ii in 2:k) XI[ii - 1, ii] = tempB

  return(XI)
}

Y_matrix_gen = function(k, nobs, y) {

  XI = XI_gen(k = k)

  Y_matrix = matrix(data = 0.0, nrow = nobs, ncol = k - 1)

  for (ii in 1:nobs) {Y_matrix[ii, ] = XI[, y[ii]]}

  return(Y_matrix)

}

beta_kernel = function(y, k, my, warm, lambda, inv_LK){

  nobs = length(y)
  dnobs = as.double(nobs)
  tnobs = NROW(inv_LK)
  beta = matrix(data = 0, nrow = tnobs, ncol = (k - 1))
  beta0 = matrix(data = 0, nrow = 1, ncol = (k - 1))

  for(q in 1:(k - 1)) {
    temp = numeric(nobs)
    temp0 = 0
    for( ii in 1L:nobs ) {
      for( jj in 1:k ) {
        t1 = warm[ii, jj] * my[jj, q]
        t2 = (2 * {y[ii] == jj} - 1)

        temp[ii] = temp[ii] + t2 * t1
        temp0 = temp0 + t2 * t1
      }
    }

    beta[, q] = inv_LK %*% (c(temp, rep(0, tnobs - nobs)))
    beta0[, q] = temp0 / dnobs / lambda
  }
  beta0 = as.vector(beta0)
  # beta0 = matrix(colSums(beta), nrow = 1)
  return(list(beta = beta, beta0 = beta0))
}


make_knn_graph_mat = function(X, k = 6)
{
  distance = fields::rdist(X, X)
  #    distance = sv.kernel(X.mat, X.mat, kernel = list(type="rbf", par=2))
  #    distance = 1/distance
  # distance[distance < 1e-6] = 0
  knn_mat = matrix(0, nrow(X), nrow(X))
  order_mat = apply(distance, 2, order)
  for(i in 1:ncol(knn_mat)) {
    knn_mat[order_mat[1:(k + 1), i], i] = 1
  }
  graph_mat = matrix(0, nrow(X), nrow(X))
  graph_mat[(t(knn_mat) + knn_mat) != 0] = 1
  diag(graph_mat) = 0
  return(graph_mat)
}

make_L_mat = function(X, kernel = "radial", kparam = 1, graph, weightType = c("Heatmap", "Binary"), normalized = FALSE)
{
  # make edge weights matrix W
  W_mat = main_kernel(X, X, kernel = list(type = kernel, par = kparam))
  W_mat = W_mat * graph
  if(weightType == "Binary")
  {
    # print("binary weight")
    W_mat = graph
  }
  d = rowSums(W_mat)
  D_mat = diag(d, nrow(W_mat), ncol(W_mat))
  L_mat = D_mat - W_mat

  if (normalized) {
    # standardize L_mat
    std_D_mat = diag(1 / sqrt(d))
    L_mat = std_D_mat %*% L_mat %*% std_D_mat
  }

  # binary matrix
  return(L_mat)
}



# predict_kernel = function(K_test, beta, beta0, k)
# {
#   n = nrow(K_test)
#
#   XI = XI_gen(k = k)
#
#   beta0 = matrix(beta0,
#                  nrow = n,
#                  ncol = ncol(beta),
#                  byrow = TRUE)
#
#   f_matrix = t(K_test %*% beta + beta0)
#
#   inner_matrix = matrix(data = 0, nrow = n, ncol = k)
#
#   for(ii in 1:k) inner_matrix[, ii] = colSums(f_matrix * XI[, ii])
#
#   z = apply(X = inner_matrix, MARGIN = 1, FUN = pred)
#
#   return(list(class = z, inner_prod = inner_matrix))
#
# }

# pred = function(f) {
#   tst = sapply(f, function(i) {isTRUE(all.equal(i, max(f)))})
#   y = min(which(tst))
#   return(y)
# }


make_anovaKernel = function(x, u, kernel)
{
  if (!is.matrix(x))  # degenerate case: x is a row vector
  { x = t(as.matrix(x))}
  else { x = as.matrix(x)}

  u = as.matrix(u)
  dimx = ncol(x)

  # calculate anova kernels for main effects
  if(kernel$type == "spline")
  {
    # assign the number of anova kernels
    numK = 2*dimx
    # list of kernel matrices
    anova_kernel = vector(mode="list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode="list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }
  }
  else if (kernel$type == 'spline2')
  {
    numK = (2*dimx) + (2*dimx*(2*dimx-1)/2 - dimx)
    anova_kernel = vector(mode="list", numK)
    kernelCoord = vector(mode="list", numK)
    index = 0
    # main effects
    for(d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }
    # two-way interactions
    for (i in 1:(dimx-1))
    {
      for (j in (i+1):dimx)
      {
        index = index + 1
        A.linear = as.matrix(anova_kernel[[2*i-1]])
        A.smooth = as.matrix(anova_kernel[[2*i]])
        B.linear = as.matrix(anova_kernel[[2*j-1]])
        B.smooth = as.matrix(anova_kernel[[2*j]])
        anova_kernel[[index]] = A.linear*B.linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep="")
        index = index + 1
        anova_kernel[[index]] = A.linear*B.smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep="")
        index = index + 1
        anova_kernel[[index]] = A.smooth*B.linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep="")
        index = index + 1
        anova_kernel[[index]] = A.smooth*B.smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep="")
      }
    }
  }
  else if (kernel$type == "spline-t")
  {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
  }
  else if (kernel$type == 'spline-t2')
  {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[, d])
      B = as.matrix(u[, d])
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
    for (i in 1:(dimx - 1))
    {
      for (j in (i + 1):dimx)
      {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep="")
      }
    }
  } else if (kernel$type == "radial2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx)
    {
      index = index + 1
      A = as.matrix(x[,d])
      B = as.matrix(u[,d])
      anova_kernel[[index]] = main_kernel(A, B, kernel)
      kernelCoord[[index]] = paste("x", d, sep="")
    }
    for (i in 1:(dimx - 1))
    {
      for (j in (i + 1):dimx)
      {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep="")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)

    # anova_kernel = mclapply(1:dimx, function(d) {main_kernel(as.matrix(x[, d]), as.matrix(u[, d]), kernel)}, mc.cores = 10)

    # for (d in 1:dimx)
    # {
    #   kernelCoord[[d]] = paste("x", d, sep = "")
    # }

     for (d in 1:dimx)
     {
       A = as.matrix(x[, d])
       B = as.matrix(u[, d])
       anova_kernel[[d]] = main_kernel(A, B, kernel)
       kernelCoord[[d]] = paste("x", d, sep = "")
     }
  }
  list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel)
}

combine_kernel = function(anova_kernel, theta = rep(1, anova_kernel$numK))
{
  K = 0
  for (d in 1:anova_kernel$numK)
  {
    K = (K + theta[d] * anova_kernel$K[[d]])
  }
  return(K)
}

class_code = function(y, k = max(y))
{
  # k: the number of classes
  classcodes = diag(rep(k / (k-1), k)) - 1 / (k-1) * matrix(1, k, k)

  # degenerate case
  if (length(y) == 1) {
    Y = t(as.matrix(classcodes[y, ]))
  } else {
    Y = as.matrix(classcodes[y, ])
  }
  return(Y)
}

find_nonzero = function(Amat)
{
  nr = nrow(Amat)
  nc = ncol(Amat)
  Amat_compact = matrix(0, nr, nc)
  Aind = matrix(0, nr + 1, nc)
  for (j in 1:nc) {
    index = (1:nr)[Amat[, j] != 0]
    number = length(index)
    Amat_compact[1:number, j] = Amat[index, j]
    Aind[1, j] = number
    Aind[2:(number+1), j] = index
  }
  max_number = max(Aind[1, ])
  Amat_compact = Amat_compact[1:max_number, ]
  Aind = Aind[1:(max_number + 1), ]
  return(list(Amat_compact = Amat_compact, Aind = Aind))
}

code = function(y)
{
  n_class = length(unique(y))
  n = length(y)
  y_index = cbind(1:n, y)

  In = diag(n)
  Lmat = matrix(1, nrow = n, ncol = n_class)
  Lmat[y_index] = 0
  Emat = diag(n_class)

  vmatj = vector("list", length = n_class)
  umatj = vector("list", length = n_class)
  for (k in 1:n_class) {
    lvecj = Lmat[, k]
    evecj = Emat[k, , drop = FALSE]
    vmatj[[k]] = diag(lvecj)
    umatj[[k]] = evecj %x% In
  }

  AH = matrix(0, n, n * n_class)
  for (k in 1:n_class) {
    AH = AH + (2 * vmatj[[k]] - In) %*% umatj[[k]] / n_class
  }

  Hmatj = vector("list", length = n_class)
  for (k in 1:n_class) {
    Hmatj[[k]] = (In - 2 * vmatj[[k]]) %*% umatj[[k]] + AH
  }

  return(list(In = In, vmatj = vmatj, umatj = umatj, AH = AH, Hmatj = Hmatj, y_index = y_index))
}

# adjacency_knn = function(X, distance = "euclidean", k = 6)
# {
#     Ds = as.matrix(dist(X, method = distance))
#     neighbours = apply(Ds, 1, function(x) sort(x, index.return = TRUE)$ix[2:(k + 1)])
# 	neighbours = as.integer(neighbours)
#     adj = as.matrix(Matrix::sparseMatrix(i = rep(1:nrow(X), each = k), j = neighbours, x = 1, dims = c(nrow(X), nrow(X))))
#     adj = (adj | t(adj)) * 1
# 	return(adj)
# }

# fixit = function(A, epsilon = .Machine$double.eps)
# {
#   eig = eigen(A, symmetric = TRUE)
#   n = length(eig$values)
#   eps = epsilon * abs(eig$values[1])
#   eig$values[eig$values < eps] = eps
#   return(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
# }

fixit = function(A, epsilon = 100 * .Machine$double.eps, is_diag = FALSE)
{
  if (is_diag) {
    d = diag(A)
    # n = length(d)
    # tol = n * epsilon
    tol = epsilon
    eps = tol * max(d)
    # if (any(d < eps)) {
    #   d = d - min(d) + eps
    # }
    d[d < eps] = eps
    Q = diag(d)
  } else {
    eig = eigen(A, symmetric = TRUE)
    n = length(eig$values)
    # tol = n * epsilon
    tol = epsilon
    eps = tol * abs(eig$values[1])
    if (any(eig$values < eps)) {
      eig$values[eig$values < eps] = eps
    }
    Q = eig$vectors %*% (eig$values * t(eig$vectors))
  }
  return(Q)
}

inverse = function(A, epsilon = 100 * .Machine$double.eps)
{
  eig = eigen(A, symmetric = TRUE)
  # n = length(eig$values)
  # tol = n * epsilon
  tol = epsilon
  eps = max(tol * abs(eig$values[1]), 0)
  positive = eig$values > eps
  if (all(positive)) {
    Q = eig$vectors %*% (1 / eig$values * t(eig$vectors))
  } else {
    Q = eig$vectors[, positive, drop = FALSE] %*% ((1 / eig$values[positive]) * t(eig$vectors[, positive, drop = FALSE]))
  }
  return(Q)
}


# fixit = function(A, epsilon)
# {
#   eig = eigen(A, symmetric = TRUE)
#   eps = epsilon * eig$values[1]
#   if (eig$values[length(eig$values)] < eps) {
#     delta = eps - eig$values[length(eig$values)]
#   } else {
#     delta = 0
#   }
#   # eig$values[eig$values < eps] = eps
#   # eig$values[eig$values <= eps] = eig$values[eig$values <= eps] + delta
#   eig$values = eig$values + delta
#   return(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
# }

# fixit = function(A, epsilon)
# {
#   d = dim(A)[1]
#   if (dim(A)[2] != d)
#     stop("Input matrix is not square!")
#   es = eigen(A, symmetric = TRUE)
#   esv = es$values
#   if (missing(epsilon)) {
#     epsilon = d * max(abs(esv)) * .Machine$double.eps
#   }
#   delta = 2 * d * max(abs(esv)) * epsilon
#   tau = pmax(0, delta - esv)
#   dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
#   return(A + dm)
# }



# fixit = function(A, epsilon)
# {
#   eig = eigen(A, symmetric = TRUE)
#   # eps = ncol(A) * epsilon * abs(eig$values[1])
#   eps = 0
#   eig$values[eig$values < eps] = eps
#   # eig$values[eig$values <= epsilon] = eig$values[eig$values <= epsilon] + epsilon
#   A = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
#   diag(A) = diag(A) + epsilon
#   return(A)
# }

# fixit = function(A, epsilon)
# {
#   return(nearPD(A)$mat)
# }
