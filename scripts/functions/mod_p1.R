mod_p1 <- function(in.data, 
                   in.vars, 
                   in.id_var, 
                   in.iter = 1000,
                   in.seed = 1234, 
                   METHOD = 0, 
                   REMOVE_OUTLIER = "Y",
                   cov.est.method = "ST", 
                   UCL_pct = 99,
                   out_boot_data = "N", 
                   in.parallel = "no",
                   weight = NULL) {
  
  # Phase I Multivariate Outlier Detection to get Upper Control Limit (UCL)
  
  # Args:
  #    in.data: Input data frame. This should not have any missing values
  #    in.vars: A character vector contains all numeric variables 
  #    in.id_var: A variable to uniquely identify an observation in the ip data. 
  #               This is the variable on the x axis of the control chart
  #    in.iter: Number of iterations for bootstrapping
  #    in.seed: A random seed used for bootstrapping
  #    METHOD: A number. 0 is for no bootstrapping, 
  #            1 is for modifed bootstrapping method for small sample size 
  #            use cases, 2 is for standard bootstrapping method using the boot 
  #            library which allows for parellel processing
  #    REMOVE_OUTLIER: 'Y' to remove outlier and 'N' to keep all
  #    cov.est.method: 'ST' for standard covaraiance estimation method (default) 
  #                    'SD' for succesive difference method
  #    UCL_pct: percentile for UCL. e.g.: 99 is for 99% UCL
  #    out_boot_data: 'Y' to output bootstrapped data and 'N' otherwise. 
  #                   Only used when METHOD is set to 1
  #    in.parallel: The type of parallel operation to be used (if any) by the 
  #                 boot function from the boot package
  #    weight: a weight vector used to weight the parametrics
  #    outliersRemove: 'Y' to remove outlier and 'N' to keep all
  #    uclPerct: Percentile for UCL. e.g.: 99 is for 99% UCL
  
  # Returns:
  #    s: Sample variance-covariance matrix used for scoring
  #    inv_S: Inverse of the sample variance-covariance matrix
  #    inv_S_i: A list of matrices which are inverse of the sample 
  #              variance-covariance matrix excluding the ith variable
  #    weight: weight matrix  
  #    xbar: Mean vector of the input data frame
  #    UCL: Upper control limit
  #    UCL_pct: Percentile for UCL
  #    model.data: the raw data for Phase I mod
  #    model.var: a vector of all the variables used for Phase I mod
  #    t2_outlier_quantile: A number indiciating percentage of removed outliers
  #    id.var: A variable to uniquely identify an observation in input data
  #    message: description for warning messages


  ######## check input parameters #########################
  message <- "The program ran successfully."

  if (!is.null(weight)) {
    if (length(in.vars) != length(weight)) {
      message <- paste("The number of weights(", length(weight), ") 
                       does not match the number of variables(",
        length(in.vars), ").",
        sep = ""
      )
      warning(message)
      return(message = message)
    }
  }

  ############### end of check input parameters #####################

  if (!data.table::is.data.table(in.data)) {
    in.data <- data.table::setDT(in.data)
  }


  if (REMOVE_OUTLIER == "Y") {

    # finding mahalanobis distance:
    xbar <- colMeans(in.data[, in.vars, with = FALSE])

    if (cov.est.method == "ST") {
      xcov <- cov(in.data[, in.vars, with = FALSE])
    } else if (cov.est.method == "SD") {
      V <- apply(in.data[, in.vars, with = FALSE], 2, diff)
      xcov <- t(V) %*% V / 2 / nrow(V)
    }

    inv_xcov <- tryCatch(solve(xcov), error = function(e) {
      return(NULL)
    })

    if (is.null(inv_xcov)) {
      inv_xcov <- tryCatch(MASS::ginv(xcov), error = function(e) {
        return(NULL)
      })
    }

    if (is.null(inv_xcov)) {
      message <- "Singular matrix. Remove some columns"
      warning(message)
      return(message = message)
    }

    if (is.null(weight)) {
      my_mahal <- mahalanobis(
        x = in.data[, in.vars, with = FALSE], center = xbar,
        cov = inv_xcov, inverted = TRUE
      )
    } else {

      # weighted mahalanobis distance
      weighted.inv_xcov <- diag(weight) %*% inv_xcov
      my_mahal <- mahalanobis(
        x = in.data[, in.vars, with = FALSE], center = xbar,
        cov = weighted.inv_xcov, inverted = TRUE
      )
    }

    in.data <- cbind(in.data, mahal = my_mahal)

    rules <- boxplot_outliers(dt = my_mahal, input_b = 4, iqr_mult = 2)
    orig_num_row <- nrow(in.data)

    in.data <- in.data[mahal < rules$adjusted_rule]
    new_num_row <- nrow(in.data)
    in.data[, mahal := NULL]
    t2_outlier_quantile <- 1 - new_num_row / orig_num_row
  } else {
    t2_outlier_quantile <- 0
  }

  # calculating xbar, S, inv_S, and inv_S_i

  xdat <- in.data[, in.vars, with = FALSE]
  out.xbar <- colMeans(xdat)

  if (cov.est.method == "ST") {
    out.S <- cov(xdat)
  } else if (cov.est.method == "SD") {
    V <- apply(xdat, 2, diff)
    out.S <- t(V) %*% V / 2 / nrow(V)
  }

  out.inv_S <- tryCatch(solve(out.S), error = function(e) {
    return(NULL)
  })

  if (is.null(out.inv_S)) {
    out.inv_S <- tryCatch(MASS::ginv(out.S), error = function(e) {
      return(NULL)
    })
  }

  if (is.null(out.inv_S)) {
    message <- "Singular matrix. Remove some columns"
    warning(message)
    return(message = message)
  }

  sample_size <- nrow(in.data)
  p <- length(in.vars)

  inv_S_i <- list()

  for (i in 1:p) {
    S_i <- out.S[-i, -i]

    if (exists("tmp.inv_S_i", inherits = FALSE)) {
      rm(tmp.inv_S_i)
    }

    tmp.inv_S_i <- tryCatch(solve(S_i), error = function(e) {
      return(NULL)
    })

    if (is.null(tmp.inv_S_i)) {
      tmp.inv_S_i <- tryCatch(MASS::ginv(S_i), error = function(e) {
        return(NULL)
      })
    }

    if (is.null(tmp.inv_S_i)) {
      message <- "Singular matrix. Remove some columns"
      warning(message)
      return(message = message)
    }

    inv_S_i[[i]] <- tmp.inv_S_i
  }

  # calculate UCL

  if (METHOD == 0) {
    deltaX <- as.matrix(xdat - out.xbar[col(xdat)])

    # calculate Hotelling T-square based on all observations and 
    # variance-covariance matrix from phase 1 data

    if (is.null(weight)) {
      T2 <- colSums(t(deltaX %*% out.inv_S) * t(deltaX))
    } else {
      weighted.inv_S <- diag(weight) %*% out.inv_S
      T2 <- colSums(t(deltaX %*% weighted.inv_S) * t(deltaX))
    }

    invisible(gc())

    UCL <- quantile(T2, UCL_pct / 100)

    return(list(
      S = out.S, inv_S = out.inv_S, inv_S_i = inv_S_i, weight = weight,
      xbar = out.xbar, UCL = UCL, UCL_pct = UCL_pct,
      model.data = xdat, model.var = in.vars,
      t2_outlier_quantile = t2_outlier_quantile,
      id.var = in.id_var, boot_data = NULL, message = message
    ))
  } else if (METHOD == 1) {
    if (!is.null(in.seed)) {
      set.seed(in.seed)
    }

    for (b in 1:in.iter) {

      # sample with replacement
      train <- sample(sample_size, size = sample_size, replace = TRUE)

      sample <- in.data[train]

      options(warn = -1) # turn warning off for duplicated samples
      oob <- in.data[-train]
      options(warn = 0)

      # calculate inverse of sample variance-covariance matrix from sample
      M <- sample[, in.vars, with = FALSE]

      if (cov.est.method == "ST") {
        sigma_hat <- cov(M)
      } else if (cov.est.method == "SD") {
        V <- apply(M, 2, diff)
        sigma_hat <- t(V) %*% V / 2 / nrow(V)
      }

      inv_sigma_hat <- tryCatch(solve(sigma_hat), error = function(e) {
        return(NULL)
      })

      if (is.null(inv_sigma_hat)) {
        inv_sigma_hat <- tryCatch(MASS::ginv(sigma_hat), error = function(e) {
          return(NULL)
        })
      }


      if (is.null(inv_sigma_hat)) {
        message <- "Singular matrix. Remove some columns"
        warning(message)
        return(message = message)
      }

      invisible(gc())

      xbar <- colMeans(M)

      # get OOB drives
      O <- oob[, in.vars, with = FALSE]
      deltaX <- as.matrix(O - xbar[col(O)])

      # calculate Hotelling T-square based on observations from oob and 
      # variance-covariance matrix from sample

      if (is.null(weight)) {
        T2 <- colSums(t(deltaX %*% inv_sigma_hat) * t(deltaX))
      } else {
        weighted.inv_sigma_hat <- diag(weight) %*% inv_sigma_hat
        T2 <- colSums(t(deltaX %*% weighted.inv_sigma_hat) * t(deltaX))
      }

      invisible(gc())

      oob[, T2 := T2][, iter := b]

      # stack all Hotelling T-squares
      if (b == 1) out.data <- data.table::copy(oob) else out.data <- 
        data.table::rbindlist(list(out.data, oob), use.names = TRUE)
    }

    UCL <- quantile(out.data$T2, UCL_pct / 100)

    if (out_boot_data == "Y") {
      out_boot_data <- out.data
    } else {
      out_boot_data <- NULL
    }

    return(list(
      S = out.S, inv_S = out.inv_S, inv_S_i = inv_S_i, weight = weight,
      xbar = out.xbar, UCL = UCL, UCL_pct = UCL_pct,
      model.data = xdat, model.var = in.vars,
      t2_outlier_quantile = t2_outlier_quantile,
      id.var = in.id_var, boot_data = out_boot_data,
      message = message
    ))
  } else if (METHOD == 2) {

    # standard bootstrapping  using library boot, allow parallel processing

    if (!is.null(in.seed)) {
      set.seed(in.seed)
    }

    test <- tryCatch(boot::boot(
      data = in.data, statistic = hotelling_t2, R = in.iter, in.vars = in.vars,
      UCL_pct = UCL_pct, cov.est.method = cov.est.method, 
      parallel = in.parallel, weight = weight
    ), error = function(e) {
      # warning(e)
      return(NULL)
    })

    if (is.null(test)) {
      message <- "Error running bootstrapping"
      warning(message)
      return(message = message)
    }

    out.T2 <- data.frame(test$t)
    names(out.T2) <- paste("P", UCL_pct, sep = "")

    UCL <- colMeans(out.T2) # use mean of percentile

    if (out_boot_data == "Y") {
      out_boot_data <- out.T2
    } else {
      out_boot_data <- NULL
    }

    return(list(
      S = out.S, inv_S = out.inv_S, inv_S_i = inv_S_i, weight = weight,
      xbar = out.xbar, UCL = UCL, UCL_pct = UCL_pct,
      model.data = xdat, model.var = in.vars,
      t2_outlier_quantile = t2_outlier_quantile,
      id.var = in.id_var, boot_data = out_boot_data, message = message
    ))
  }
}
