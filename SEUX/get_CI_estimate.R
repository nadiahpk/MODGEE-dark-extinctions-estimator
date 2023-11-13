get_CI_estimate <- function (S, E, nreps = 10000, percentile = 95, U_T = 0, midP_fnc = NULL) 
  {
    if (!is.vector(S)) {
        stop("Input S needs to be a vector.")
    }
    if (!is.vector(E)) {
        stop("Input E needs to be a vector.")
    }
    if (!is.numeric(U_T)) {
        stop("Input U_T needs to be a whole number.")
    }
    else {
        if (U_T%%1 != 0) {
            stop("Input U_T needs to be a whole number: ", as.character(U_T))
        }
    }
    if (is.null(midP_fnc)) {
        midP_fnc <- central_hyper_midP
    }
    else {
        if (typeof(midP_fnc) != "closure") {
            stop("midP_fnc needs to a function.")
        }
    }
    if (nreps%%1 != 0) {
        stop("Input nreps needs to be whole no: ", as.character(nreps))
    }
    if (!is.numeric(percentile)) {
        stop("Percentile needs to be a percentage number: ", 
            as.character(percentile))
    }
    if (percentile < 1) {
        stop("Percentile needs to be a percentage (e.g. 95): ", 
            as.character(percentile))
    }
    if (length(E) != length(S)) {
        stop("Input vectors S and E must be have length.")
    }
    T_idx <- length(S)
    d <- S[-1] - S[-T_idx] + E[-1] - E[-T_idx]
    U_min_possible <- c(rev(cumsum(rev(d[1:length(d)]))), 0)
    UM <- matrix(NA, nreps, T_idx)
    XM <- matrix(NA, nreps, T_idx)
    for (nrep in 1:nreps) {
        U <- replicate(T_idx, 0)
        U[T_idx] <- U_T
        impossibleFlag <- FALSE
        for (t in T_idx:2) {
            alpha <- runif(1)
            S0 <- S[t - 1]
            S1 <- S[t]
            d0 <- d[t - 1]
            U1 <- U[t]
            res <- find_U0_bnd(midP_fnc, alpha, S0, S1, U1, d0, 
                impossibleFlag)
            U[t - 1] <- res$U0
            impossibleFlag <- res$impossibleFlag
        }
        U <- pmax(U, U_min_possible)
        N <- S[1] + E[1] + U[1]
        X <- N - E - S - U
        X <- pmax(X, 0)
        UM[nrep, ] <- U
        XM[nrep, ] <- X
    }
    q_lo <- (1 - percentile/100)/2
    q_hi <- 1 - (1 - percentile/100)/2
    X_lo <- apply(XM, 2, function(x) quantile(x, q_lo))
    X_hi <- apply(XM, 2, function(x) quantile(x, q_hi))
    X_mean <- apply(XM, 2, function(x) mean(x))
    U_lo <- apply(UM, 2, function(x) quantile(x, q_lo))
    U_hi <- apply(UM, 2, function(x) quantile(x, q_hi))
    U_mean <- apply(UM, 2, function(x) mean(x))
    return(list(data.frame(X_lo = X_lo, X_hi = X_hi, X_mean = X_mean, 
        U_lo = U_lo, U_hi = U_hi, U_mean = U_mean),XM=XM,UM=UM))
  }