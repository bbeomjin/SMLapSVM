smlapsvm = function(x = NULL, y, gamma = 0.5, ux = NULL, valid_x = NULL, valid_y = NULL, nfolds = 5, type = c("rm", "ram"),
                    lambda_seq = 2^{seq(-10, 10, length.out = 100)}, lambda_I_seq = 2^{seq(-20, 15, length.out = 20)},
                    lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                    adjacency_k = 6, normalized = FALSE, weightType = "Binary",
                    kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                    scale = FALSE, criterion = c("0-1", "loss"), isCombined = TRUE, optModel = FALSE, nCores = 1, verbose = 0, ...)
{
  out = list()
  call = match.call()
  type = match.arg(type)
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  if (type == "ram") {
    cstep_fun = cstep.sramlapsvm
    thetastep_fun = thetastep.sramlapsvm
  } else {
    if (gamma == 0) {
      cstep_fun = cstep.smlapsvm
      thetastep_fun = thetastep.smlapsvm
    } else {
      cstep_fun = cstep.srmlapsvm
      thetastep_fun = thetastep.srmlapsvm
    }
  }

  cat("Fit c-step \n")
  cstep_args = list(x = x, y = y, gamma = gamma, ux = ux, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, lambda_I_seq = lambda_I_seq,
                    theta = NULL, adjacency_k = adjacency_k, normalized = normalized, weightType = weightType,
                    kernel = kernel, kparam = kparam, scale = scale, criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  if ((gamma == 0) & (type == "rm")) {
    cstep_args$gamma = NULL
  }

  cstep_fit = do.call(cstep_fun, cstep_args)

  cat("Fit theta-step \n")

  if (optModel) {
    thetastep_opt = FALSE
  } else {
    thetastep_opt = TRUE
  }

  thetastep_fit = thetastep_fun(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined, nCores = nCores, ...)

  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
    cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
  }

  out$opt_theta = thetastep_fit$opt_theta
  out$cstep_inform = cstep_fit
  out$thetastep_inform = thetastep_fit

  if (optModel) {
    cat("Fit c-step \n")
    cstep_args$theta = thetastep_fit$opt_theta
    cstep_args$optModel = TRUE
    opt_cstep_fit = do.call(cstep_fun, cstep_args)

    if (verbose == 1) {
      cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
    }

    out$opt_model = opt_cstep_fit$opt_model
    out$opt_valid_err = opt_cstep_fit$opt_valid_err
    out$valid_err = opt_cstep_fit$valid_err
  }
  out$type = type
  out$call = call
  class(out) = "smlapsvm"
  return(out)
}

predict.smlapsvm = function(object, newx = NULL)
{
  if (is.null(newx)) {
    newx = object$cstep_inform$x
  }

  if (is.null(object$opt_model)) {
    model = object$thetastep_inform$opt_model
  } else {
    model = object$opt_model
  }

  new_anova_K = make_anovaKernel(newx, rbind(object$cstep_inform$x, object$cstep_inform$ux),
                                 kernel = object$cstep_inform$kernel, object$cstep_inform$kparam)
  newK = combine_kernel(new_anova_K, object$opt_theta)

  if (object$type == "rm") {
    pred = predict.rmlapsvm_compact(model, newK = newK)
  } else {
    pred = predict.ramlapsvm_compact(model, newK = newK)
  }
  return(list(class = pred$class, pred_value = pred$pred_value))
}



smsvm = function(x = NULL, y, gamma = 0.5, valid_x = NULL, valid_y = NULL, nfolds = 5, type = c("rm", "ram"),
                 lambda_seq = 2^{seq(-10, 10, length.out = 100)},
                 lambda_theta_seq = 2^{seq(-10, 10, length.out = 100)},
                 kernel = c("linear", "gaussian", "poly", "spline", "anova_gaussian"), kparam = c(1),
                 scale = FALSE, criterion = c("0-1", "loss"),
                 isCombined = TRUE, optModel = FALSE, nCores = 1, verbose = 0, ...)
{
  out = list()
  call = match.call()
  type = match.arg(type)
  kernel = match.arg(kernel)
  criterion = match.arg(criterion)

  if (type == "ram") {
    cstep_fun = cstep.sramsvm
    thetastep_fun = thetastep.sramsvm
  } else {
    if (gamma == 0) {
      cstep_fun = cstep.smsvm
      thetastep_fun = thetastep.smsvm
    } else {
      cstep_fun = cstep.srmsvm
      thetastep_fun = thetastep.srmsvm
    }
  }

  cat("Fit c-step \n")
  cstep_args = list(x = x, y = y, gamma = gamma, valid_x = valid_x, valid_y = valid_y, nfolds = nfolds,
                    lambda_seq = lambda_seq, theta = NULL, kernel = kernel, kparam = kparam,
                    criterion = criterion, optModel = FALSE, nCores = nCores, ...)

  if ((gamma == 0) & (type == "rm")) {
    cstep_args$gamma = NULL
  }

  cstep_fit = do.call(cstep_fun, cstep_args)

  cat("Fit theta-step \n")

  if (optModel) {
    thetastep_opt = FALSE
  } else {
    thetastep_opt = TRUE
  }

  thetastep_fit = thetastep_fun(cstep_fit, lambda_theta_seq = lambda_theta_seq, isCombined = isCombined,
                                optModel = thetastep_opt, nCores = nCores, ...)

  if (verbose == 1) {
    cat("CV-error(cstep):", cstep_fit$opt_valid_err, "\n")
    cat("CV-error(theta-step):", thetastep_fit$opt_valid_err, "\n")
  }

  out$opt_theta = thetastep_fit$opt_theta
  out$cstep_inform = cstep_fit
  out$thetastep_inform = thetastep_fit

  if (optModel) {
    cat("Fit c-step \n")
    cstep_args$theta = thetastep_fit$opt_theta
    cstep_args$optModel = TRUE
    opt_cstep_fit = do.call(cstep_fun, cstep_args)

    if (verbose == 1) {
      cat("CV-error(cstep):", opt_cstep_fit$opt_valid_err, "\n")
    }

    out$opt_model = opt_cstep_fit$opt_model
    out$opt_valid_err = opt_cstep_fit$opt_valid_err
    out$valid_err = opt_cstep_fit$valid_err
  }
  out$type = type
  out$call = call
  class(out) = "smsvm"
  return(out)
}


predict.smsvm = function(object, newx = NULL)
{
  if (is.null(newx)) {
    newx = object$cstep_inform$x
  }

  if (is.null(object$opt_model)) {
    model = object$thetastep_inform$opt_model
  } else {
    model = object$opt_model
  }

  new_anova_K = make_anovaKernel(newx, object$cstep_inform$x,
                                 kernel = object$cstep_inform$kernel, object$cstep_inform$kparam)
  newK = combine_kernel(new_anova_K, object$opt_theta)

  if (object$type == "rm") {
    pred = predict.rmsvm_compact(model, newK = newK)
  } else {
    pred = predict.ramsvm_compact(model, newK = newK)
  }
  return(list(class = pred$class, pred_value = pred$pred_value))
}


