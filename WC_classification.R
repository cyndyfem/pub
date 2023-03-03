supersom<-function (data, grid = somgrid(), rlen = 100, alpha = c(0.05, 
    0.01), radius = quantile(nhbrdist, 2/3), whatmap = NULL, 
    user.weights = 1, maxNA.fraction = 0L, keep.data = TRUE, 
    dist.fcts = NULL, mode = c("online", "batch", "pbatch"), 
    cores = -1, init, normalizeDataLayers = TRUE) 
{
    data <- check.data(data)
    narows <- check.data.na(data, maxNA.fraction = maxNA.fraction)
    data <- remove.data.na(data, narows)
    full.data <- data
    nmat <- length(data)
    whatmap <- check.whatmap(data, whatmap)
    nmap <- length(whatmap)
    data <- data[whatmap]
    grid <- check.somgrid(grid)
    nhbrdist <- unit.distances(grid)
    if (length(radius) == 1) 
        radius <- c(radius, 0)
    nobjects <- nrow(data[[1]])
    nvar <- sapply(data, ncol)
    data.matrix <- matrix(unlist(data), ncol = nobjects, byrow = TRUE)
    nNA <- getnNA(data, maxNA.fraction, nobjects)
    if (length(dist.fcts) == length(full.data)) {
        orig.dist.fcts <- dist.fcts
        dist.fcts <- dist.fcts[whatmap]
    }
    else {
        if (length(dist.fcts) == 1) {
            orig.dist.fcts <- rep(dist.fcts, length(full.data))
            dist.fcts <- orig.dist.fcts[whatmap]
        }
        else {
            defaultDist <- "sumofsquares"
            defaultClassDist <- "tanimoto"
            factorLayers <- sapply(full.data, function(datamat) is.factor(datamat) || 
                is.factor.matrix(datamat))
            default.dist.fcts <- rep(defaultDist, length(full.data))
            if (any(factorLayers)) 
                default.dist.fcts[which(factorLayers)] <- defaultClassDist
            if (length(dist.fcts) == 0) {
                orig.dist.fcts <- default.dist.fcts
                dist.fcts <- orig.dist.fcts[whatmap]
            }
            else {
                if (length(dist.fcts) == length(whatmap)) {
                  orig.dist.fcts <- default.dist.fcts
                  orig.dist.fcts[whatmap] <- dist.fcts
                }
                else {
                  stop("Wrong number of distances defined")
                }
            }
        }
    }
    if (any(tanidists <- dist.fcts == "tanimoto")) {
        minvals <- sapply(data[which(tanidists)], min, na.rm = TRUE)
        maxvals <- sapply(data[which(tanidists)], max, na.rm = TRUE)
        if (any(minvals < 0 | maxvals > 1)) 
            stop("Layers for which the Tanimoto distance is used should have data within the [0,1] range")
    }
    dist.ptrs <- getDistancePointers(dist.fcts, maxNA.fraction = maxNA.fraction)
    ncodes <- nrow(grid$pts)
    if (missing(init)) {
        starters <- sample(1:nobjects, ncodes, replace = FALSE)
        init <- lapply(data, function(x) x[starters, , drop = FALSE])
    }
    else {
        if (!is.list(init) & (is.vector(init) | is.factor(init) | 
            is.matrix(init))) 
            init <- list(init)
        if (length(init) != nmap) 
            stop("Incorrect number of initialization matrices")
        if (!all(sapply(init, nrow) == ncodes)) 
            stop("Incorrect number of objects in initalization matrices")
        if (!all(sapply(init, ncol) == nvar)) {
            if (maxNA.fraction == 0L) {
                stop("Incorrect number of variables in initialization matrices")
            }
            else {
                stop("Incorrect number of variables in initialization matrices, ", 
                  "maybe due to the removal of columns because of NAs")
            }
        }
    }
    if (maxNA.fraction > 0L) {
        for (jj in seq(along = init)) {
            nastarters <- which(is.na(init[[jj]]), arr.ind = TRUE)
            if (nrow(nastarters) > 0) {
                for (i in unique(nastarters[, 1])) {
                  idx <- which(nastarters[, 1] == i)
                  imputers <- sample(data[[jj]][i, !is.na(data[[jj]][i, 
                    ])], length(idx), replace = TRUE)
                  init[[jj]][i, nastarters[idx, 2]] <- imputers
                }
            }
        }
    }
    init.matrix <- matrix(unlist(init), ncol = ncodes, byrow = TRUE)
    distance.weights <- orig.user.weights <- rep(0, nmat)
    if (length(whatmap) == 1) {
        weights <- user.weights <- 1
        distance.weights[whatmap] <- orig.user.weights[whatmap] <- 1
    }
    else {
        if (length(user.weights) == 1) {
            user.weights <- rep(user.weights, length(whatmap))
        }
        else {
            if (length(user.weights) == nmat) 
                user.weights <- user.weights[whatmap]
        }
        if (any(user.weights == 0)) 
            stop("Incompatibility between whatmap and user.weights")
        if (abs(sum(user.weights)) < .Machine$double.eps) 
            stop("user.weights sum to zero")
        user.weights <- user.weights/sum(user.weights)
        orig.user.weights[whatmap] <- user.weights
        if (normalizeDataLayers) {
            meanDistances <- lapply(seq(along = init), function(ii) object.distances(list(data = init[ii], 
                whatmap = 1, user.weights = 1, distance.weights = 1, 
                maxNA.fraction = maxNA.fraction, dist.fcts = dist.fcts[ii]), 
                type = "data"))
            if (any(sapply(meanDistances, mean) < .Machine$double.eps)) 
                stop("Non-informative layers present: mean distance between objects zero")
            distance.weights[whatmap] <- 1/sapply(meanDistances, 
                mean)
        }
        else {
            distance.weights <- rep(1, length(data))
        }
        weights <- user.weights * distance.weights[whatmap]
        weights <- weights/sum(weights)
    }
    mode <- match.arg(mode)
    switch(mode, online = {
        res <- RcppSupersom(data = data.matrix, codes = init.matrix, 
            numVars = nvar, weights = weights, distanceFunctions = dist.ptrs, 
            numNAs = nNA, neighbourhoodDistances = nhbrdist, 
            neighbourhoodFct = as.integer(grid$neighbourhood.fct), 
            alphas = alpha, radii = radius, numEpochs = rlen)
    }, batch = {
        res <- RcppBatchSupersom(data = data.matrix, codes = init.matrix, 
            numVars = nvar, weights = weights, distanceFunctions = dist.ptrs, 
            numNAs = nNA, neighbourhoodDistances = nhbrdist, 
            neighbourhoodFct = as.integer(grid$neighbourhood.fct), 
            radii = radius, numEpochs = rlen)
    }, pbatch = {
        res <- RcppParallelBatchSupersom(data = data.matrix, 
            codes = init.matrix, numVars = nvar, weights = weights, 
            distanceFunctions = dist.ptrs, numNAs = nNA, neighbourhoodDistances = nhbrdist, 
            neighbourhoodFct = as.integer(grid$neighbourhood.fct), 
            radii = radius, numEpochs = rlen, numCores = cores)
    })
    changes <- matrix(res$changes, ncol = nmap, byrow = TRUE)
    colnames(changes) <- names(data)
    mycodes <- res$codes
    layerID <- rep(1:nmap, nvar)
    mycodes2 <- split(as.data.frame(mycodes), layerID)
    mycodes3 <- lapply(mycodes2, function(x) t(as.matrix(x)))
    codes <- vector(length(full.data), mode = "list")
    names(codes) <- names(full.data)
    codes[whatmap] <- mycodes3
    for (ii in seq(along = whatmap)) colnames(codes[[whatmap[ii]]]) <- colnames(data[[ii]])
    if (keep.data) {
        mapping <- map.kohonen(structure(list(codes = codes, 
            distance.weights = distance.weights, dist.fcts = orig.dist.fcts, 
            data = full.data), class = "kohonen"), whatmap = whatmap, 
            user.weights = orig.user.weights, maxNA.fraction = maxNA.fraction)
        structure(list(data = full.data, unit.classif = mapping$unit.classif, 
            distances = mapping$distances, grid = grid, codes = codes, 
            changes = changes, alpha = alpha, radius = radius, 
            na.rows = narows, user.weights = orig.user.weights, 
            distance.weights = distance.weights, whatmap = whatmap, 
            maxNA.fraction = maxNA.fraction, dist.fcts = orig.dist.fcts), 
            class = "kohonen")
    }
    else {
        structure(list(grid = grid, codes = codes, changes = changes, 
            alpha = alpha, radius = radius, na.rows = narows, 
            user.weights = orig.user.weights, distance.weights = distance.weights, 
            whatmap = whatmap, maxNA.fraction = maxNA.fraction, 
            dist.fcts = orig.dist.fcts), class = "kohonen")
    }
}

yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"),
                           whatmap = c("alpha", "cdc15", "cdc28", "elu"),
                           maxNA.fraction = .5)


somgrid <-function (xdim = 8, ydim = 6, topo = c("rectangular", "hexagonal"), 
    neighbourhood.fct = c("bubble", "gaussian"), toroidal = FALSE) 
{
    topo <- match.arg(topo)
    x <- 1L:xdim
    y <- 1L:ydim
    pts <- as.matrix(expand.grid(x = x, y = y))
    if (topo == "hexagonal") {
        pts[, 1L] <- pts[, 1L] + 0.5 * (pts[, 2L]%%2)
        pts[, 2L] <- sqrt(3)/2 * pts[, 2L]
    }
    neighbourhood.fct <- match.arg(neighbourhood.fct)
    neighbourhood.fct <- factor(neighbourhood.fct, levels = c("bubble", 
        "gaussian"))
    res <- list(pts = pts, xdim = xdim, ydim = ydim, topo = topo, 
        neighbourhood.fct = neighbourhood.fct, toroidal = toroidal)
    class(res) <- "somgrid"
    res
}

mygrid <- somgrid(7, 7, "hexagonal")
WSesom <- list(grid = mygrid)
class(WSsom) <- "kohonen"

WS.supersom <- supersom(WS, somgrid(7, 7, "hexagonal"),
                           whatmap = c("alpha", "cdc15", "cdc28", "elu"),
                           maxNA.fraction = .5)
