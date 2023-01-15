sigest = function(x, q)
{
  a = quantile(as.numeric(dist(x, method = "euclidean")), q)
  return(1 / (2 * a^2))
}

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
  } else if(kernel == "gaussian" | kernel == "gaussian-2way") {
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

spline_kernel = function(x, u)
{
  x <- as.matrix(x)
  u <- as.matrix(u)
  K1x <- (x - 1/2)
  K1u <- (u - 1/2)
  K2x <- (K1x^2 - 1/12)/2
  K2u <- (K1u^2 - 1/12)/2
  ax <- x%x%matrix(1, 1, nrow(u))
  au <- u%x%matrix(1, 1, nrow(x))
  b <- abs(ax - t(au))
  K1 <- K1x%x%t(K1u)
  K2 <- K2x%x%t(K2u) - ((b - 1/2)^4 - (b - 1/2)^2/2 + 7/240)/24
  list(K1 = K1, K2 = K2)
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
  } else if (kernel == "gaussian-2way") {
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

code_rmsvm2 = function(y)
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



fixit = function(A, epsilon = .Machine$double.eps) {

  if (!is.matrix(A)) {
    A = as.matrix(A)
  }
  d = dim(A)
  eig = eigen(A, symmetric = TRUE)
  v = eig$values
  delta = max(abs(v)) * epsilon
  tau = pmax(0, delta - v)
  A = eig$vectors %*% diag(v + tau, d[1]) %*% t(eig$vectors)
  # A = eig$vectors %*% ((v + tau) * t(eig$vectors))
  return(A)
}

prediction_err = function(y, pred, type = "0-1")
{
  if (type == "0-1") {
    # standard accuracy
    acc = sum(y == pred) / length(y)
  } else if (type == "balanced") {
    # balanced accuracy
    y = factor(y, levels = unique(c(y, pred)))
    pred = factor(pred, levels = levels(y))
    t = table(y, pred)
    temp = diag(t) / rowSums(t)
    temp[is.nan(temp)] = 0
    acc = mean(temp)
  } else {
    # not yet
  }
  return(acc)
}





