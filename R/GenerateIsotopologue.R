#' @export
setMethod(
  'GenerateIsotopologue',
  signature = c(param = 'IsotopologueParam', experimentParam = 'ExperimentParam'),
  function(param = IsotopologueParam(),
           experimentParam = ExperimentParam()) {
    wd0 <- getwd()
    setwd(experimentParam@wd)
    obj <- new('IsotopologueTargets')
    obj@.processHistory <- append(obj@.processHistory,
                                  list('GenerateIsotopologue' = param))
    obj@ExperimentParam <- experimentParam

    if (!file.exists(param@file.id)) {
      stop(paste0("File '", param@file.id, "' could not be found in '",
                  experimentParam@wd, "'!"))
    }
    if (!file.exists(param@file.xset)) {
      stop(paste0("File '", param@file.xset, "' could not be found in '",
                  experimentParam@wd, "'!"))
    }

    xset <- LoadData(param@file.xset)

    dt.anno.raw <- read.csv(param@file.id, stringsAsFactors = FALSE)
    dt.anno <- dt.anno.raw[which(!is.na(dt.anno.raw$ID)), , drop = FALSE]
    rownames(dt.anno) <- dt.anno$name
    col.info <- c('mz', 'isotope', 'adduct', 'compound.name', 'ID', 'Formula')
    info.list <- lapply(dt.anno$name, function(nm) {
      dr <- dt.anno[nm, ]
      info <- do.call(cbind, strsplit(as.character(dr[col.info]), split = ';'))
      colnames(info) <- col.info
      info <- info[!grepl('\\+\\d+', info[, 'isotope']), , drop = FALSE]
      if (nrow(info) > 0) {
        return(info)
      } else {
        return(NULL)
      }
    })
    names(info.list) <- dt.anno$name
    idx.keep <- which(!sapply(info.list, is.null))
    info.list <- info.list[idx.keep]
    dt.anno <- dt.anno[idx.keep, , drop = FALSE]
    isotope.list.pre <- BiocParallel::bplapply(info.list, function(info) {
      require(MetTracer)
      res <- lapply(seq(nrow(info)), function(nr) {
        dr <- info[nr, ]
        element.trace <- sapply(gsub('\\d+', '', experimentParam@element.trace),
                                GetElementNumber, dr['Formula'])
        if (sum(element.trace) == 0) {
          return(NA)
        }
        res.iso <- GetIsotopePattern(formula = dr['Formula'],
                                     adduct = dr['adduct'],
                                     element.trace = experimentParam@element.trace,
                                     element.label = experimentParam@element.label,
                                     ppm = experimentParam@ppm,
                                     res.define = experimentParam@res.define,
                                     is.plot.pattern = param@is.plot.pattern)
        list('info' = c('formula' = res.iso$formula,
                        dr[c('compound.name', 'ID', 'adduct')]),
             'iso' = res.iso$iso)
      })
      names(res) <- info[, 'ID']
      if (all(is.na(res))) {
        return(NA)
      }
      res <- res[!is.na(res)]
      iso <- lapply(res, `[[`, 'iso')

      info.sum <- t(sapply(res, `[[`, 'info'))
      return(list('iso' = iso, 'info' = info.sum))
    })
    idx.keep <- which(!is.na(isotope.list.pre))
    isotope.list.pre <- isotope.list.pre[idx.keep]
    dt.anno <- dt.anno[idx.keep, , drop = FALSE]

    isotope.list <- lapply(isotope.list.pre, `[[`, 'iso')
    info.formula <- lapply(isotope.list.pre, `[[`, 'info')

    groupidx <- match(names(isotope.list), xcms::groupnames(xset))
    eic.present <- GetFeatureEIC(xset, groupidx, extra = param@rt.extend,
                                 method = param@method, value = param@value)
    mz.isomer <- lapply(isotope.list, function(isotope) {
      sort(unique(round(do.call(rbind, isotope)[, 1], 6)))
    })
    rt.rg.eic <- t(sapply(eic.present, function(eic) {
      range(eic[, 1])
    }))
    colnames(rt.rg.eic) <- c('rtmin', 'rtmax')
    mz.rg.eic <- lapply(mz.isomer, function(isotope) {
      mz <- sort(unique(unlist(isotope)))
      t(sapply(mz, PpmRange, experimentParam@ppm))
    })

    smp.present <- attributes(eic.present)$smp.present
    dt.anno$smp.present <- xcms::sampnames(xset)[smp.present]
    dt.anno$maxo.present <- GetPresentValue(xset, "maxo", smp.present, groupidx)
    dt.anno$into.present <- GetPresentValue(xset, "into", smp.present, groupidx)
    dt.anno$intb.present <- GetPresentValue(xset, "intb", smp.present, groupidx)
    dt.anno$top3.present <- sapply(eic.present, function(eic) {
      idx.apex <- which(eic[, "inpeak"] == 2)
      sum(eic[seq(idx.apex - 1, idx.apex + 1), "intensity"])
    })

    dt.value <- xcms::groupval(xset, value = param@value)
    rownames(dt.value) <- xcms::groupnames(xset)
    dt.anno <- cbind(dt.anno, dt.value[rownames(dt.anno), , drop = FALSE])

    obj@TargetInfo <- dt.anno
    obj@Isotopologue <- isotope.list
    obj@Formula <- info.formula
    obj@Profile <- eic.present
    obj@rtRange <- rt.rg.eic
    obj@mzRange <- mz.rg.eic
    obj@mzIsomer <- mz.isomer
    write.csv(dt.anno, file = file.path(experimentParam@d.tmp, "PeakInfo.csv"),
              row.names = FALSE)
    saveRDS(obj, file = file.path(experimentParam@d.tmp, 'IsotopologueTargets.rda'),
            version = 2)
    setwd(wd0)

    return(obj)
  })

setGeneric('GetIsotopePattern', function(formula, adduct,
                                         element.trace = '12C',
                                         element.label = '13C',
                                         ppm = 25,
                                         res.define = 400,
                                         is.plot.pattern = FALSE) {

  name.element.label <- gsub('^\\d+', '', element.label)
  num.element.extra <- rep(0, length(element.label))
  names(name.element.label) <- names(num.element.extra) <- element.trace

  data(adducts, package = getPackageName())
  data(elements, package = getPackageName())
  data(isotopes, package = 'enviPat')
  formula.checked <- enviPat::check_chemform(isotopes, formula)[,2]

  adduct.info <- adducts[which(adducts$Name == adduct), ]
  if (nrow(adduct.info) == 0) stop(paste0("Adduct type '", adduct, "' not found!"))
  formula.adduct <- enviPat::multiform(formula.checked, adduct.info$Multi)

  if (adduct.info$Formula_add != 'FALSE') {
    num.element.extra <- sapply(name.element.label, GetElementNumber,
                                adduct.info$Formula_add)
    formula.adduct <- enviPat::mergeform(formula.adduct, adduct.info$Formula_add)
  }

  if (adduct.info$Formula_ded != 'FALSE') {
    formula.adduct <- enviPat::subform(formula.adduct, adduct.info$Formula_ded)
  }

  elements.labeled <- rbind(elements,
                            isotopes[which(isotopes$isotope %in%
                                             c(element.label, element.trace)), ,
                                     drop = FALSE])
  elements.labeled <- unique(elements.labeled[order(elements.labeled$element), , drop = FALSE])

  pattern <- enviPat::isopattern(elements.labeled,
                                 formula.adduct,
                                 threshold = 0,
                                 plotit = is.plot.pattern,
                                 charge = adduct.info$Charge,
                                 emass = 0.00054858,
                                 verbose = FALSE)[[1]]

  is.avaliable <- sapply(element.trace, function(element) {
    pattern[, element] >= num.element.extra[element]
  })
  row.keep <- which(apply(is.avaliable, 1, all))
  if (length(element.label) == 1) {
    profiles <- enviPat::envelope(list(pattern[row.keep, ]), verbose = FALSE)
  } else {
    mz.mono <- min(pattern[row.keep, 1])
    dmz = max(mz.mono, res.define) * ppm * 1e-6
    profiles <- enviPat::envelope(list(pattern[row.keep, ]), dmz = dmz,
                                  verbose = FALSE)
  }
  centro <- enviPat::vdetect(profiles, plotit = is.plot.pattern,
                             verbose = FALSE)
  return(list('iso' = centro[[1]],
              'formula' = formula.adduct))
})

setGeneric('GetFeatureEIC', function(xset, groupidx, extra = 30,
                                     method = 'max', value = 'into') {
  if (missing(groupidx)) {
    is.all <- TRUE
  }
  if (!all(file.exists(xset@filepaths))) {
    files <- xset@filepaths
    basename(files)
    fp.all <- list.files(pattern = '.mzXML$', recursive = TRUE, full.names = TRUE)
    xset@filepaths <- fp.all[match(basename(files), basename(fp.all))]
  }

  groupname <- xcms::groupnames(xset)

  rt.range.data <- range(xset@rt)

  rt.min.origin <- xcms::groupval(xset, value = 'rtmin')
  rt.max.origin <- xcms::groupval(xset, value = 'rtmax')
  rt.apex.origin <- xcms::groupval(xset, value = 'rt')
  rownames(rt.apex.origin) <- rownames(rt.min.origin) <- rownames(rt.max.origin) <- groupname
  rtmin <- apply(rt.min.origin, 1, median) - extra
  rtmax <- apply(rt.max.origin, 1, median) + extra

  rtmin[which(rtmin < rt.range.data[1])] <- rt.range.data[1]
  rtmax[which(rtmax > rt.range.data[2])] <- rt.range.data[2]

  rt.range.eic <- cbind(rtmin, rtmax)

  if (missing(groupidx)) {
    groupidx <- seq(nrow(rt.range.eic))
  }

  eics <- BiocParallel::bplapply(seq(length(xset@filepaths)), function(idx) {
    xcms::getEIC(xset, rtrange = rt.range.eic[groupidx, , drop = FALSE],
                 sampleidx = idx,
                 groupidx = groupidx)
  })
  eics.profile <- sapply(eics, function(eic) {
    names(eic@eic[[1]]) <- eic@groupnames
    eic@eic[1]
  })
  smpname <- names(eics.profile)

  eic.present <- switch(
    method,
    'max' = {
      idx.max <- apply(xcms::groupval(xset, value = value), 1, which.max)
      names(idx.max) <- groupname
      rt.apex <- sapply(groupname[groupidx], function(nm) {
        rt.apex.origin[nm, idx.max[nm]]
      })

      res <- sapply(eics[[1]]@groupnames, function(nm) {
        eic <- eics.profile[[idx.max[nm]]][[nm]]
        idx.s <- sum(eic[, 'rt'] < rt.min.origin[nm, idx.max[nm]])
        idx.e <- sum(eic[, 'rt'] <= rt.max.origin[nm, idx.max[nm]])
        idx.apex <- which.min(abs(eic[, 'rt'] - rt.apex[nm]))
        inpeak <- eic[, 'rt'] >= rt.min.origin[nm, idx.max[nm]] &
          eic[, 'rt'] <= rt.max.origin[nm, idx.max[nm]]
        inpeak[idx.apex] <- 2
        cbind(eic, inpeak)
      }, simplify = FALSE)
      attr(res, "smp.present") <- idx.max[groupidx]
      return(res)
    })
  return(eic.present)
})

setGeneric(
  "GetElementNumber",
  function(element, formula) {
    reg.expr <- paste0('(?<=', element, ')\\d+')
    num.element <- regmatches(formula,
                              gregexpr(reg.expr, formula,
                                       perl = TRUE))[[1]]
    ifelse(length(num.element), as.numeric(num.element), 0)
  })

setGeneric(
  "GetPresentValue",
  function(xset, val.name, smp.present, groupidx) {
    dt.val <- xcms::groupval(xset, value = val.name)[groupidx, , drop = FALSE]
    val <- apply(cbind(seq(length(groupidx)), smp.present), 1, function(dr) {
      dt.val[dr[1], dr[2]]
    })
    return(val)
  })

