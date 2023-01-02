# make a matrix symmetric
makeSymm <- function(m, from_lower = T) {
    if (from_lower)
        m[upper.tri(m)] = t(m)[upper.tri(m)] # from lower
    else
        m[lower.tri(m)] = t(m)[lower.tri(m)] # from upper
    return(m)
}