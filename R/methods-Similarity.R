setGeneric(
  "GetSimilarity",
  function(peak.target, peak.extract, rt) {
    apex.extract <- which.min(abs(peak.extract$rt - rt))
    apex.target <- which(peak.target$inpeak == 2)
    lside <- min(apex.extract, apex.target) - 1
    rside <- min(nrow(peak.extract) - apex.extract,
                 nrow(peak.target) - apex.target)
    int.target <- peak.target$intensity.s[seq(from = apex.target - lside,
                                              to   = apex.target + rside)]
    int.extract <- peak.extract$intensity.s[seq(from = apex.extract - lside,
                                                to   = apex.extract + rside)]
    res <- cor(int.target, int.extract, method="pearson")
    res <- ifelse(is.na(res), 0, res)
    return(res)
  })
