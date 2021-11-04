#' @export
setMethod(
  'ExtractIsotopicPeaks',
  signature = c(param.extract = 'ExtractParam',
                param.pd = 'PeakdetectionParam',
                iso.targets = 'IsotopologueTargets'),
  function(param.extract = ExtractParam(),
           param.pd = PeakdetectionParam(),
           iso.targets,
           cutoff = 0.6,
           check.peaks = TRUE,
           ...) {

    wd0 <- getwd()

    iso.peaks <- new("IsotopicPeaks")
    iso.peaks@.processHistory <- list(iso.targets@.processHistory,
                                      "ExtractIsotopicPeaks" =
                                        list("param.extract" = param.extract,
                                             "param.pd" = param.pd))
    iso.peaks@ExperimentParam <- experimentParam <- iso.targets@ExperimentParam
    iso.peaks@TargetInfo <- iso.targets@TargetInfo

    setwd(experimentParam@wd)
    files <- list.files(param.extract@d.extract, pattern = '(?i).mzxml$',
                        recursive = TRUE, full.names = TRUE)
    smp.groups <- xcms:::phenoDataFromPaths(files)
    if (ncol(smp.groups) > 1) {
      smp.groups$class <- apply(smp.groups, 1, paste, collapse = "/")
      smp.groups <- smp.groups[, -(seq(ncol(smp.groups) - 1)), drop = FALSE]
    } else {
      smp.groups$class <- as.character(smp.groups$class)
    }
    smp.groups$files <- files
    smp.groups$samples <- seq_along(files)
    smp.groups$samplenames <- rownames(smp.groups)
    iso.peaks@SampleGroups <- smp.groups

    mz.isomer <- iso.targets@mzIsomer
    mzrange <- iso.targets@mzRange
    rtrange <- iso.targets@rtRange

    require(ggsci)
    cls <- pal_d3("category20")(20)[-c(3,4)]
    ncl <- length(cls)

    nms.ft <- names(iso.targets@Isotopologue)
    names(nms.ft) <- nms.ft
    method.best <- param.extract@method.best

    ####################################
    data.extract <- BiocParallel::bplapply(files, function(file) {
      require(MetTracer)
      d.plot <- file.path(experimentParam@d.tmp, 'PeakCheck', method.best)
      if (!dir.exists(d.plot)) {
        dir.create(d.plot, recursive = TRUE)
      }
      if (check.peaks) {
        pdf(file.path(d.plot, paste0(basename(file), ".pdf")),
            height = 12, width = 8)
        par(mfrow = c(3,1))
      }

      xr <- xcms::xcmsRaw(file)

      peaks.extract <- switch(
        param.extract@method.align,
        'apex' = {
          eics <- ExtractIsotopeEIC(xr, mzrange, mz.isomer)
          param.mspk <- msPeaks:::ParamFindPeaks(methodRoi = param.pd@method.roi,
                                                 methodBaseline = param.pd@method.baseline,
                                                 methodSmooth = param.pd@method.smooth)
          param <- switch(param.pd@method.peakdetection,
                          'centWave' = {
                            msPeaks::ParamCentWave(peakwidth = param.pd@peakwidth,
                                                   snthr = param.pd@snthr,
                                                   prefilter = param.pd@prefilter,
                                                   fitgauss = param.pd@fitgauss,
                                                   ParamFindPeaks = param.mspk)
                          },
                          'localMax' = {
                            msPeaks::ParamLocalMax(peakwidth = param.pd@peakwidth,
                                                   snthr = param.pd@snthr,
                                                   smoothNoise = param.pd@smoothNoise,
                                                   ParamFindPeaks = param.mspk)
                          })


          peaks.ft <- lapply(nms.ft, function(nm.ft) {
            eic.target <- iso.targets@Profile[[nm.ft]]
            eic.target <- msPeaks::SmoothEIC(eic = data.frame(eic.target),
                                             peakwidth = param.pd@peakwidth,
                                             method = param.pd@method.smooth)
            eic.target$intensity.s[is.na(eic.target$intensity.s)] <- 0
            peak.target <- eic.target[eic.target[, "inpeak"] > 0, ]
            eics.mz <- eics[[nm.ft]]
            peaks.mz <- lapply(names(eics.mz), function(nm.mz) {
              eic <- eics.mz[[nm.mz]]
              peaks <- msPeaks::FindPeaks(param,
                                          intensity = eic$intensity,
                                          rt = eic$rt,
                                          rtrange = rtrange[nm.ft, ])
              if (nrow(peaks@peaks) == 0) {
                return(NA)
              }
              rownames(peaks@peaks) <- paste(nm.mz, rownames(peaks@peaks),
                                             sep = '_')
              names(peaks@eics) <- rownames(peaks@peaks)
              peaks
            })

            names(peaks.mz) <- names(eics.mz)

            peaks.roi <- lapply(peaks.mz, function(peaks) {
              if (class(peaks) == "MSPeaks") {
                peaks@roi
              } else {
                NA
              }
            })
            is.find.roi <- !is.na(peaks.roi)
            peaks.roi <- peaks.roi[is.find.roi]
            peaks.mz <- peaks.mz[!is.na(peaks.mz)]

            if (length(peaks.mz) == 0) {
              return(NA)
            }
            peaks.table <- lapply(names(peaks.mz), function(nm) {
              peaks.mz[[nm]]@peaks
            })
            peaks.table <- do.call(rbind, peaks.table)
            peaks.eics <- lapply(peaks.mz, function(peaks) {
              peaks@eics
            })
            names(peaks.eics) <- NULL
            peaks.eics <- do.call(c, peaks.eics)

            cor.pearson <- sapply(rownames(peaks.table), function(nm.pk) {
              GetSimilarity(peak.target,
                            peaks.eics[[nm.pk]],
                            peaks.table[nm.pk, "rt"])
            })

            peaks.table$cor <- cor.pearson
            is.ok <- cor.pearson > cutoff
            peaks.table$best <- 0

            if (any(is.ok)) {
              nm.best <- do.call(paste0('FindBestPeak.', method.best),
                                 list("peaks.table" = peaks.table,
                                      "peak.target" = peak.target,
                                      "cutoff" = cutoff,
                                      "check.cluster" = check.peaks))
              peaks.table[nm.best, "best"] <- 1

              nm.best.eic <- strsplit(nm.best, split = '_')[[1]][1]
              rt.iso <- peaks.table[nm.best, 'rt']

              eic.best <- eics.mz[[nm.best.eic]]
              peak.best <- peaks.eics[[nm.best]]
              apex.eic.best <- which.min(abs(eic.best$rt - rt.iso))
              apex.peak.best <- which.min(abs(peak.best$rt - rt.iso))

              lside <- apex.peak.best - 1
              rside <- nrow(peak.best) - apex.peak.best
              rg.iso <- seq(from = apex.eic.best - lside,
                            to   = apex.eic.best + rside)


              eic.iso <- lapply(eics.mz, function(eic) {
                eic <- eic[rg.iso, ]
              })
              names(eic.iso) <- names(eics.mz)
              baseline.iso <- sapply(names(eic.iso), function(nm.eic) {
                roi <- peaks.roi[[nm.eic]]
                if (is.null(roi)) {
                  return(NA)
                }
                nr <- which(roi[, "scmin"] <= apex.eic.best &
                              roi[, "scmax"] >= apex.eic.best)
                if (length(nr) == 0) {
                  return(NA)
                }
                return(roi[nr, 'baseline'])
              })
            } else {
              return(NA)
            }

            if (check.peaks) {
              if (any(!is.na(eic.iso))) {

                plot(eic.target[, 1:2], ylim = c(0, max(c(eic.target[, 2],
                                                          peaks.table$maxo))),
                     type = "b", col = "gray", main = nm.ft)
                lines(peak.target[, c(1, 4)], col = 'green4', lty = 2, lwd = 2)
                lapply(seq_along(peaks.mz), function(idx) {
                  lapply(peaks.mz[[idx]]@eics, function(eic) {
                    lines(eic[, c(1, 3)], col = cls[idx %% ncl])
                  })
                })
                legend("topleft", names(peaks.mz),
                       col = cls[seq(length(peaks.mz))],
                       lty = 1,
                       cex = 0.8)


                plot(eic.target[, 1:2], ylim = c(0, max(c(eic.target[, 2],
                                                          peaks.table$maxo))),
                     type = "b", col = "gray", main = nm.ft)
                lines(peak.target[, c(1, 4)], col = 'green4', lty = 2, lwd = 2)
                lwds <- (names(eic.iso) %in% nm.best.eic*2 + 1)
                lapply(seq_along(eic.iso), function(idx) {
                  lines(eic.iso[[idx]], col = cls[idx %% ncl],
                        lwd = lwds[idx])
                })
                legend("topleft", names(eic.iso),
                       col = cls[seq(length(eic.iso))],
                       lty = 1,
                       lwd = lwds)
              }
            }
            return(list("peaks.iso" = peaks.table,
                        "eics.iso" = eic.iso,
                        "info.iso" = data.frame('baseline' = baseline.iso,
                                                'findroi' = is.find.roi),
                        "apex.iso" = apex.peak.best))
          })
          peaks.ft
        },
        'eic' = {
          eics <- ExtractIsotopeROI(xr, mzrange, mz.isomer)
        }
      )

      if (check.peaks) {
        dev.off()
      }
      return(peaks.extract)
    })

    names(data.extract) <- rownames(smp.groups)
    saveRDS(data.extract,
            file = file.path(experimentParam@d.tmp,
                             paste0("isoExtract_", method.best, ".rda")))

    if (check.peaks) {
      d.sv <- file.path(experimentParam@d.res, method.best)
      dir.create(d.sv, recursive = T)

      lapply(names(data.extract), function(nm.smp) {
        data.smp <- data.extract[[nm.smp]]
        data.smp <- data.smp[!is.na(data.smp)]
        pt <- lapply(names(data.smp),function(nm.ft) {
          pk <- data.smp[[nm.ft]]
          pk.iso <- pk$peaks.iso
          pk.iso$nm <- nm.ft
          pk.iso
        })
        pt <- do.call(rbind, pt)

        write.csv(pt, file = file.path(d.sv, paste0(basename(nm.smp), '.csv')))
      })
    }

    iso.peaks@Peaks <- plyr::ldply(nms.ft, function(nm.ft) {
      do.call(rbind, lapply(seq_along(data.extract), function(idx.smp) {
        data.smp <- data.extract[[idx.smp]][[nm.ft]]
        if (any(!is.null(data.smp)) & all(!is.na(data.smp))) {
          pks <- data.smp$peaks.iso
          pks$samples <- idx.smp
          pks$peakname <- rownames(pks)
          pks$featurename <- nm.ft
        } else {
          pks <- NULL
        }
        pks
      }))
    })
    iso.profile <- lapply(nms.ft, function(nm.ft) {
      nms.mz <- as.character(mz.isomer[[nm.ft]])
      profile.mz <- lapply(nms.mz, function(nm.mz) {
        profile.smp <- lapply(seq_along(data.extract), function(idx.smp) {
          data.smp <- data.extract[[idx.smp]][[nm.ft]]
          if (any(!is.null(data.smp)) & all(!is.na(data.smp))) {
            data.smp[['eics.iso']][[nm.mz]]
          } else {
            NA
          }
        })
        names(profile.smp) <- smp.groups$samplenames
        profile.smp
      })
      names(profile.mz) <- nms.mz
      profile.mz
    })

    iso.apex <- t(sapply(nms.ft, function(nm.ft) {
      sapply(seq_along(data.extract), function(idx.smp) {
        data.smp <- data.extract[[idx.smp]][[nm.ft]]
        if (any(!is.null(data.smp)) & all(!is.na(data.smp))) {
          apex <- data.smp$apex.iso
        } else {
          apex <- NA
        }
        apex
      })
    }))

    colnames(iso.apex) <- smp.groups$samplenames

    iso.info <- lapply(nms.ft, function(nm.ft) {
      info <- lapply(seq_along(data.extract), function(idx.smp) {
        data.smp <- data.extract[[idx.smp]][[nm.ft]]
        if (any(!is.null(data.smp)) & all(!is.na(data.smp))) {
          data.smp$info.iso
        } else {
          NA
        }
      })
      names(info) <- smp.groups$samplenames
      info
    })

    iso.peaks@Profile <- iso.profile
    iso.peaks@Apex <- data.frame(iso.apex)
    iso.peaks@Info <- iso.info

    saveRDS(iso.peaks,
            file = file.path(experimentParam@d.tmp,
                             paste0("isoPeaks_", method.best, ".rda")))
    setwd(wd0)

    return(iso.peaks)
  })

#' Extract EIC of isomers
#'
#' Extract EIC of isomers according to their mz range cross all scans and
#' return the mzs and intensities of the EIC as well as the scantime for the
#' file
#'
#' @param xr `xcmsSet` xcmsSet dataset from xcms
#' @param mzrange `list` mz ranges for extraction
#' @param mz.isomer `list` mz values for isomers
#' @param rtrange `list` rt ranges for extraction
#'
#' @return list EICs
#' @export
setGeneric(
  'ExtractIsotopeEIC',
  function(xr, mzrange, mz.isomer = NULL, rtrange = NULL) {
    rtrange.scan <- range(xr@scantime)
    eics <- sapply(names(mzrange), function(nm.ft) {
      mzrange.ft <- mzrange[[nm.ft]]
      if (is.null(rtrange)) {
        rtrange.ft <- rtrange.scan
      } else {
        rtrange.ft <- rtrange[nm.ft, ]
      }
      eics.ft <- plyr::alply(mzrange.ft, 1, function(r.mzrange.ft) {
        eic <- xcms::rawEIC(xr, r.mzrange.ft, rtrange.ft)
        eic <- data.frame('rt' = xr@scantime[eic$scan],
                          'intensity' = eic$intensity)
      })

      if (!is.null(mz.isomer)) {
        names(eics.ft) <- mz.isomer[[nm.ft]]
      }
      return(eics.ft)
    })
    return(eics)
  })

setGeneric(
  'ExtractIsotopeROI',
  function(file, iso.targets, ppm = 25) {
    xr <- xcms::xcmsRaw(file)
    rt.rg.raw <- range(xr@scantime)
    info.targets <- iso.targets@TargetInfo
    rownames(info.targets) <- info.targets$name
    isotopelogue <- iso.targets@Isotopologue
    profile <- iso.targets@Profile
    rt.rg <- iso.targets@rtRange
    mz.rg <- iso.targets@mzRange

    eics.iso <- sapply(names(mz.rg), function(nm) {
      rt.apex <- profile[[nm]][which(profile[[nm]][, 'inpeak'] == 2), 'rt']
      idx.apex <- which.min(abs(xr@scantime - rt.apex))
      mz.rg.1 <- mz.rg[[nm]]
      Roi <- apply(mz.rg.1, 1, function(mz.rg.ion) {
        eic <- xcms::rawEIC(xr, mz.rg.ion, rt.rg.raw)
        is.roi <- getContinuousPtsAboveThrIdx(eic$intensity, 0, 4,
                                              0, 0)
        idx.fr.roi <- GetRoi(is.roi, idx.apex)
        if (!is.null(idx.fr.roi)) {
          return(do.call(cbind, eic)[idx.fr.roi, , drop = FALSE])
        } else {
          return(NULL)
        }
      })
      names(eics) <- names(isotopelogue[[nm]])
      eics
    })
  })

setGeneric(
  "FindBestPeak.maxcor",
  function(peaks.table, peak.target, check.cluster = FALSE, ...) {
    names.best <- FindBestCluster(peaks.table, peak.target, cutoff = cutoff,
                                  check.cluster = check.cluster)
    names.best[which.max(peaks.table[names.best, "cor"])]
  }
)

setGeneric(
  "FindBestPeak.maxint",
  function(peaks.table, peak.target, cutoff = 0.6, check.cluster = FALSE) {
    names.best <- FindBestCluster(peaks.table, peak.target, cutoff = cutoff,
                                  check.cluster = check.cluster)
    names.best[which.max(peaks.table[names.best, "maxo"])]
  })

setGeneric(
  "FindBestCluster",
  function(peaks.table, peak.target, cutoff = 0.6, check.cluster = FALSE) {
    return(NULL)
  }
)
