kernelMatrix = function(x, y, kernel = "gaussian", kparam = 1.0) {

  x = as.matrix(x)
  y = as.matrix(y)
  p = ncol(x)

  if (NCOL(x) == 0) {
    x = matrix(0, nrow = nrow(x), ncol = 1)
  }

  if (NCOL(y) == 0) {
    y = matrix(0, nrow = nrow(y), ncol = 1)
  }

  if (kernel == "poly") {
    K = (x %*% t(y) + 1.0)^kparam
  } else if(kernel == "gaussian" | kernel == "gaussian2") {
    normx = rowSums(x^2)
    normy = rowSums(y^2)
    temp = x %*% t(y)
    temp = (-2.0 * temp) + outer(normx, rep(1.0, nrow(y)), "*") + outer(rep(1.0, nrow(x)), normy, "*")
    K = exp(-temp * kparam)
    # K = kernlab:::kernelMatrix(rbfdot(sigma = kparam), x, y)
  } else if (kernel == "spline") {
    K = 0
    for (d in 1:p) {
      K_temp = spline_kernel(x[, d, drop = FALSE], y[, d, drop = FALSE])
      K = K + K_temp$K1 + K_temp$K2
    }
  } else if (kernel == "linear") {
    K = tcrossprod(x, y)
  } else if (kernel == "anova_gaussian") {
    K = 0
    for (d in 1:p) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = kernelMatrix(A, B, kernel = "gaussian", kparam = kparam)
      K = K + K_temp
    }
  } else {
    K = NULL
  }
  return(K)
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


# make_knn_graph_mat = function(X, k = 6)
# {
#   distance = as.matrix(dist(X, method = "euclidean"))
#   #    distance = sv.kernel(X.mat, X.mat, kernel = list(type="rbf", par=2))
#   #    distance = 1/distance
#   # distance[distance < 1e-6] = 0
#   knn_mat = matrix(0, nrow(X), nrow(X))
#   order_mat = apply(distance, 2, order)
#   for(i in 1:ncol(knn_mat)) {
#     knn_mat[order_mat[1:k, i], i] = 1
#   }
#   graph_mat = matrix(0, nrow(X), nrow(X))
#   graph_mat[(t(knn_mat) + knn_mat) != 0] = 1
#   # diag(graph_mat) = 0
#   return(graph_mat)
# }

# equal to adjacency_knn function in RSSL package
make_knn_graph_mat = function(X, distance = "euclidean", k = 6)
{
  Ds <- as.matrix(dist(X, method = distance))
  neighbours <- as.integer(apply(Ds, 1, function(x) sort(x, index.return = TRUE)$ix[2:(k + 1)]))
  adj <- as.matrix(sparseMatrix(i = rep(1:nrow(X), each = k), j = neighbours, x = 1, dims = c(nrow(X), nrow(X))))
  adj <- (adj | t(adj)) * 1
  return(adj)
}

make_L_mat = function(X, kernel = "gaussian", kparam = 1, graph, weightType = c("Heatmap", "Binary"), normalized = FALSE)
{
  # make edge weights matrix W
  W_mat = kernelMatrix(X, X, kernel = kernel, kparam = kparam)
  W_mat = W_mat * graph
  if(weightType == "Binary")
  {
    # print("binary weight")
    W_mat = graph
  }
  d = rowSums(W_mat)
  D_mat = diag(d)
  L_mat = D_mat - W_mat

  if (normalized) {
    # standardize L_mat
    std_D_mat = diag(1 / sqrt(d))
    L_mat = std_D_mat %*% L_mat %*% std_D_mat
  }

  # binary matrix
  return(L_mat)
}



make_anovaKernel = function(x, y, kernel, kparam)
{
  x = as.matrix(x)
  y = as.matrix(y)
  dimx = ncol(x)

  # calculate anova kernels for main effects
  if (kernel == "spline") {
    # assign the number of anova kernels
    numK = 2 * dimx
    # list of kernel matrices
    anova_kernel = vector(mode = "list", numK)
    # list of kernel coordinate indices
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep="")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep="")
    }

  } else if (kernel == 'spline2') {
    numK = (2 * dimx) + (2 * dimx * (2 * dimx - 1) / 2 - dimx)
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    # main effects
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = K_temp$K1
      kernelCoord[[index]] = paste("x", d, " linear", sep = "")
      index = index + 1
      anova_kernel[[index]] = K_temp$K2
      kernelCoord[[index]] = paste("x", d, " smooth", sep = "")
    }
    # two-way interactions
    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A_linear = as.matrix(anova_kernel[[2 * i - 1]])
        A_smooth = as.matrix(anova_kernel[[2 * i]])
        B_linear = as.matrix(anova_kernel[[2 * j - 1]])
        B_smooth = as.matrix(anova_kernel[[2 * j]])
        anova_kernel[[index]] = A_linear * B_linear
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_linear * B_smooth
        kernelCoord[[index]] = paste("x", i, " linear,", " x", j, " smooth", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_linear
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " linear", sep = "")
        index = index + 1
        anova_kernel[[index]] = A_smooth * B_smooth
        kernelCoord[[index]] = paste("x", i, " smooth,", " x", j, " smooth", sep = "")
      }
    }
  } else if (kernel == "spline-t") {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0
    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }
  } else if (kernel == 'spline-t2') {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      K_temp = spline_kernel(A, B)
      anova_kernel[[index]] = (K_temp$K1 + K_temp$K2)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }

    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else if (kernel == "gaussian2") {
    numK = dimx + dimx * (dimx - 1) / 2
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    index = 0

    for (d in 1:dimx) {
      index = index + 1
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[index]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[index]] = paste("x", d, sep = "")
    }

    for (i in 1:(dimx - 1)) {
      for (j in (i + 1):dimx) {
        index = index + 1
        A = anova_kernel[[i]]
        B = anova_kernel[[j]]
        anova_kernel[[index]] = A * B
        kernelCoord[[index]] = paste("x", i, " x", j, sep = "")
      }
    }
  } else {
    numK = dimx
    anova_kernel = vector(mode = "list", numK)
    kernelCoord = vector(mode = "list", numK)
    for (d in 1:dimx) {
      A = x[, d, drop = FALSE]
      B = y[, d, drop = FALSE]
      anova_kernel[[d]] = kernelMatrix(A, B, kernel, kparam)
      kernelCoord[[d]] = paste("x", d, sep = "")
    }
  }
  return(list(K = anova_kernel, coord = kernelCoord, numK = numK, kernel = kernel, kparam = kparam))
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

code_rmsvm = function(y)
{
  n_class = length(unique(y))
  n = length(y)
  y_index = cbind(1:n, y)

  In = diag(n)
  Lmat = matrix(-1, nrow = n, ncol = n_class)
  Lmat[y_index] = 1
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
    AH = AH + -vmatj[[k]] %*% umatj[[k]] / n_class
  }

  Hmatj = vector("list", length = n_class)
  for (k in 1:n_class) {
    Hmatj[[k]] = vmatj[[k]] %*% umatj[[k]] + AH
  }

  return(list(In = In, vmatj = vmatj, umatj = umatj, AH = AH, Hmatj = Hmatj, y_index = y_index))
}

code_ramsvm = function(y)
{
  n_class = length(unique(y))
  n = length(y)
  yyi = Y_matrix_gen(k = n_class, nobs = n, y = y)
  W = XI_gen(n_class)

  y_index = cbind(1:n, y)
  index_mat = matrix(-1, nrow = n, ncol = n_class)
  index_mat[y_index] = 1

  Hmatj = list()
  Lmatj = list()
  for (j in 1:(n_class - 1)) {
    Hmatj_temp = NULL
    Lmatj_temp = NULL
    for (i in 1:n_class) {
      temp = diag(n) * W[j, i]
      diag(temp) = diag(temp) * index_mat[, i]
      Hmatj_temp = rbind(Hmatj_temp, temp)
      Lmatj_temp = c(Lmatj_temp, diag(temp))
    }
    Hmatj[[j]] = Hmatj_temp
    Lmatj[[j]] = Lmatj_temp
  }
  return(list(yyi = yyi, W = W, Hmatj = Hmatj, Lmatj = Lmatj, y_index = y_index))
}

data_split = function(y, nfolds, seed = length(y))
{
  # k: the number of classes
  y = as.factor(y)
  n_data = length(y)
  n_class = length(levels(y))
  class_size = table(y)
  classname = names(class_size)
  ran = rep(0, n_data)
  if ((min(class_size) < nfolds) & (nfolds != n_data))
  {
    warning('The given fold is bigger than the smallest class size. \n Only a fold size smaller than the minimum class size \n or the same as the sample size (LOOCV) is supported.\n')
    return(NULL)
  }

  if (min(class_size) >= nfolds) {
    set.seed(seed)
    for (j in 1:n_class) {
      ran[y == classname[j]] = ceiling(sample(class_size[j]) / (class_size[j] + 1) * nfolds)
    }
  }
  else if (nfolds == n_data) {
    ran = 1:n_data
  }
  return(ran)
}



# fixit = function(A, epsilon = .Machine$double.eps) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#   d = dim(A)
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#   tol = max(abs(v)) * epsilon
#   v[v < tol] = tol
#   A = eig$vectors %*% diag(v, d[1]) %*% t(eig$vectors)
#   # A = eig$vectors %*% (v * t(eig$vectors))
#   return(A)
# }

fixit = function(A, epsilon = .Machine$double.eps) {

  if (!is.matrix(A)) {
    A = as.matrix(A)
  }

  d = dim(A)
  eig = eigen(A, symmetric = TRUE)
  # eig = eigen(A)
  v = eig$values

  # delta = 2 * d[1] * max(abs(v)) * epsilon
  delta = max(abs(v)) * epsilon

  tau = max(0, delta - v)
  A = eig$vectors %*% diag(v + tau, d[1]) %*% t(eig$vectors)
  # A = eig$vectors %*% ((v + tau) * t(eig$vectors))
  return(A)
}

# fixit4 = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
# {
#   if (is_diag) {
#     d = diag(A)
#     tol = epsilon
#     eps = max(tol * max(d), 0)
#     d[d < eps] = eps
#     A = diag(d)
#   } else {
#     eig = eigen(A, symmetric = TRUE)
#     d = eig$values
#     tol = epsilon
#     eps = max(tol * abs(d[1]), 0)
#     d[d < eps] = eps
#     # d[d < eps] = 0
#     # d = d + eps
#     Q = eig$vectors
#     o_diag = diag(A)
#     A = Q %*% (d * t(Q))
#     D = sqrt(pmax(eps, o_diag) / diag(A))
#     A[] = D * A * rep(D, each = ncol(A))
#   }
#   return(A)
# }


# fixit5 = function(A, epsilon = .Machine$double.eps) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   n = dim(A)[1]
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#   tol = max(abs(v)) * epsilon
#   diag(A) = diag(A) + tol
#   return(A)
# }


# fixit2 = function(A, epsilon = .Machine$double.eps) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   n = dim(A)[1]
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#   # tol = n * max(abs(v)) * epsilon
#   tol = max(abs(v)) * epsilon
#   positive = v > tol
#   A = eig$vectors[, positive, drop = FALSE] %*% diag(v[positive], sum(positive)) %*% t(eig$vectors[, positive, drop = FALSE])
#   diag(A) = diag(A) + tol
#   return(A)
# }


# fixit2 = function(A, epsilon = .Machine$double.eps) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   d = dim(A)
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#
#   # delta = 2 * d[1] * max(abs(v)) * epsilon
#   delta = max(abs(v)) * epsilon
#
#   tau = max(0, delta - v)
#   o_diag = diag(A)
#
#   A = eig$vectors %*% diag(v + tau, d[1]) %*% t(eig$vectors)
#   D = sqrt(o_diag / diag(A))
#   A[] = D * A * rep(D, each = ncol(A))
#   return(A)
# }

# fixit6 = function(A, epsilon = .Machine$double.eps, symm = FALSE) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   d = dim(A)
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#   tol = max(abs(v)) * epsilon
#   # tau = pmax(0, tol - v)
#   tau = pmax(0, tol - v)
#   if (symm) {
#     A = eig$vectors %*% diag(v + tau, d[1]) %*% t(eig$vectors)
#   } else {
#     eps_mat = eig$vectors %*% diag(tau, d[1]) %*% t(eig$vectors)
#     A = A + eps_mat
#   }
#   return(A)
# }

# fixit7 = function(A, epsilon = .Machine$double.eps, symm = FALSE) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   d = dim(A)
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#   tol = max(abs(v)) * epsilon
#   # tau = pmax(0, tol - v)
#   tau = pmax(0, tol - v)
#   A = eig$vectors %*% diag(v + tau, d[1]) %*% t(eig$vectors)
#   return(A)
# }

# fixit8 = function(A, epsilon = .Machine$double.eps) {
#
#   if (!is.matrix(A)) {
#     A = as.matrix(A)
#   }
#
#   d = dim(A)
#   eig = eigen(A, symmetric = TRUE)
#   # eig = eigen(A)
#   v = eig$values
#
#   delta = max(abs(v)) * epsilon
#
#   tau = pmax(0, delta + v)
#   A = eig$vectors %*% diag(tau, d[1]) %*% t(eig$vectors)
#   return(A)
# }


##############################################################################################################################################






inverse = function(A, epsilon = .Machine$double.eps)
{
  eig = eigen(A, symmetric = TRUE)
  # n = length(eig$values)
  # tol = n * epsilon
  tol = epsilon
  eps = max(tol * abs(eig$values[1]), 0)
  positive = eig$values > eps
  eig$values[!positive] = eps
  # eig$values = eig$values + eps
  Q = eig$vectors %*% ((1 / eig$values) * t(eig$vectors))
  return(Q)
}


# inverse = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
# {
#   eig = eigen(A, symmetric = TRUE)
#   tol = epsilon
#   eps = max(tol * eig$values[1], 0)
#   positive = eig$values > eps
#   Q = eig$vectors[, positive, drop = FALSE] %*% diag(1 / eig$values[positive]) %*% t(eig$vectors[, positive, drop = FALSE])
#   return(Q)
# }

# inverse = function(A, epsilon = .Machine$double.eps, is_diag = FALSE)
# {
#   d = dim(A)
#   svds = La.svd(A, d[1], d[1])
#   # tol = epsilon
#   tol = d[1] * epsilon
#   eps = tol * max(svds$d)
#   positive = svds$d > eps
#   Q = svds$u[, positive, drop = FALSE] %*% diag(1 / svds$d[positive], sum(positive)) %*% svds$vt[positive, , drop = FALSE]
#   return(Q)
# }


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

# fixit = function(A, epsilon = .Machine$double.eps)
# {
#   d = dim(A)[1]
#   if (dim(A)[2] != d)
#     stop("Input matrix is not square!")
#   es = eigen(A, symmetric = TRUE)
#   esv = es$values
#   # if (missing(epsilon)) {
#     epsilon = d * max(abs(esv)) * epsilon
#   # }
#   delta = 2 * d * max(abs(esv)) * epsilon
#   tau = pmax(0, delta - esv)
#   dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
#   A = es$vectors %*% diag(esv, d) %*% t(es$vectors)
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
