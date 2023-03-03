BC<-function (o.c, m.c, m.p, iter = 30, ratio.seq = rep(FALSE, ncol(o.c)), 
    trace = 0.05, trace.calc = 0.5 * trace, jitter.factor = 0, 
    n.tau = NULL, ratio.max = 2, ratio.max.trace = 10 * trace, 
    ties = "first", qmap.precalc = FALSE, rot.seq = NULL, silent = FALSE, 
    n.escore = 0, return.all = FALSE, subsample = NULL, pp.type = 7) 
{
    if (!is.null(rot.seq)) {
        if (length(rot.seq) != iter) {
            stop("length(rot.seq) != iter")
        }
    }
    if (length(trace.calc) == 1) 
        trace.calc <- rep(trace.calc, ncol(o.c))
    if (length(trace) == 1) 
        trace <- rep(trace, ncol(o.c))
    if (length(jitter.factor) == 1) 
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if (length(ratio.max) == 1) 
        ratio.max <- rep(ratio.max, ncol(o.c))
    if (length(ratio.max.trace) == 1) 
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    escore.iter <- rep(NA, iter + 2)
    if (n.escore > 0) {
        n.escore <- min(nrow(o.c), nrow(m.c), n.escore)
        escore.cases.o.c <- unique(suppressWarnings(matrix(seq(nrow(o.c)), 
            ncol = n.escore)[1, ]))
        escore.cases.m.c <- unique(suppressWarnings(matrix(seq(nrow(m.c)), 
            ncol = n.escore)[1, ]))
        escore.iter[1] <- escore(x = o.c[escore.cases.o.c, , 
            drop = FALSE], y = m.c[escore.cases.m.c, , drop = FALSE], 
            scale.x = TRUE)
        if (!silent) 
            cat("RAW", escore.iter[1], ": ")
    }
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if (!qmap.precalc) {
        for (i in seq(ncol(o.c))) {
            fit.qmap <- QDM(o.c = o.c[, i], m.c = m.c[, i], m.p = m.p[, 
                i], ratio = ratio.seq[i], trace.calc = trace.calc[i], 
                trace = trace[i], jitter.factor = jitter.factor[i], 
                n.tau = n.tau, ratio.max = ratio.max[i], ratio.max.trace = ratio.max.trace[i], 
                subsample = subsample, pp.type = pp.type)
            m.c.qmap[, i] <- fit.qmap$mhat.c
            m.p.qmap[, i] <- fit.qmap$mhat.p
        }
    }
    m.c <- m.c.qmap
    m.p <- m.p.qmap
    if (n.escore > 0) {
        escore.iter[2] <- escore(x = o.c[escore.cases.o.c, , 
            drop = FALSE], y = m.c[escore.cases.m.c, , drop = FALSE], 
            scale.x = TRUE)
        if (!silent) 
            cat("QDM", escore.iter[2], ": ")
    }
    m.iter <- vector("list", iter)
    o.c.mean <- colMeans(o.c)
    o.c.sdev <- apply(o.c, 2, sd)
    o.c.sdev[o.c.sdev < .Machine$double.eps] <- 1
    o.c <- scale(o.c, center = o.c.mean, scale = o.c.sdev)
    m.c.p <- rbind(m.c, m.p)
    m.c.p.mean <- colMeans(m.c.p)
    m.c.p.sdev <- apply(m.c.p, 2, sd)
    m.c.p.sdev[m.c.p.sdev < .Machine$double.eps] <- 1
    m.c.p <- scale(m.c.p, center = m.c.p.mean, scale = m.c.p.sdev)
    Xt <- rbind(o.c, m.c.p)
    for (i in seq(iter)) {
        if (!silent) 
            cat(i, "")
        if (is.null(rot.seq)) {
            rot <- rot.random(ncol(o.c))
        }
        else {
            rot <- rot.seq[[i]]
        }
        Z <- Xt %*% rot
        Z.o.c <- Z[1:nrow(o.c), , drop = FALSE]
        Z.m.c <- Z[(nrow(o.c) + 1):(nrow(o.c) + nrow(m.c)), , 
            drop = FALSE]
        Z.m.p <- Z[(nrow(o.c) + nrow(m.c) + 1):nrow(Z), , drop = FALSE]
        for (j in seq(ncol(Z))) {
            Z.qdm <- QDM(o.c = Z.o.c[, j], m.c = Z.m.c[, j], 
                m.p = Z.m.p[, j], ratio = FALSE, jitter.factor = jitter.factor[j], 
                n.tau = n.tau, pp.type = pp.type)
            Z.m.c[, j] <- Z.qdm$mhat.c
            Z.m.p[, j] <- Z.qdm$mhat.p
        }
        m.c <- Z.m.c %*% t(rot)
        m.p <- Z.m.p %*% t(rot)
        Xt <- rbind(o.c, m.c, m.p)
        if (n.escore > 0) {
            escore.iter[i + 2] <- escore(x = o.c[escore.cases.o.c, 
                , drop = FALSE], y = m.c[escore.cases.m.c, , 
                drop = FALSE], scale.x = TRUE)
            if (!silent) 
                cat(escore.iter[i + 2], ": ")
        }
        if (return.all) {
            m.c.i <- sweep(sweep(m.c, 2, attr(m.c.p, "scaled:scale"), 
                "*"), 2, attr(m.c.p, "scaled:center"), "+")
            m.p.i <- sweep(sweep(m.p, 2, attr(m.c.p, "scaled:scale"), 
                "*"), 2, attr(m.c.p, "scaled:center"), "+")
            m.iter[[i]] <- list(m.c = m.c.i, m.p = m.p.i)
        }
    }
    if (!silent) 
        cat("\n")
    m.c <- sweep(sweep(m.c, 2, m.c.p.sdev, "*"), 2, m.c.p.mean, 
        "+")
    m.p <- sweep(sweep(m.p, 2, m.c.p.sdev, "*"), 2, m.c.p.mean, 
        "+")
    for (i in seq(ncol(o.c))) {
        m.c[, i] <- sort(m.c.qmap[, i])[rank(m.c[, i], ties.method = ties)]
        m.p[, i] <- sort(m.p.qmap[, i])[rank(m.p[, i], ties.method = ties)]
    }
    names(escore.iter)[1:2] <- c("RAW", "QM")
    names(escore.iter)[-c(1:2)] <- seq(iter)
    list(mhat.c = m.c, mhat.p = m.p, escore.iter = escore.iter, 
        m.iter = m.iter)
}

out <- BC (o.c, m.c, m.p)