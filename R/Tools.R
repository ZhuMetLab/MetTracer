#' Load data to envirement by assigning it to a new object
setGeneric('LoadData', function(file, keep.name = FALSE, env){
  if (missing(env)) env <- new.env()
  b <- load(file, envir = env)
  if (keep.name | length( b) > 1) {
    r <- lapply( b, function(b1) env[[ b1]])
    names( r) <- b
    r
  } else {
    env[[b]]
  }
})

#' output a data range with given reference data and ppm
setGeneric('PpmRange', function(ref, ppm, res.def = 400) {
  dev <- ppm * 1e-6
  ref + c(-1, 1) * max(ref, res.def) * dev
})
