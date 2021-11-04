#' @export
setMethod(
  "ExtractIsotopologue",
  signature = c(param.extract = 'ExtractParam',
                iso.peaks = 'IsotopicPeaks',
                iso.targets = 'IsotopologueTargets'),
  function(param.extract = ExtractParam(),
           iso.peaks,
           iso.targets,
           ...) {
    smp.groups <- iso.peaks@SampleGroups
    peaks <- iso.peaks@Peaks
    peaks.profile <- iso.peaks@Profile
    peaks.apex <- iso.peaks@Apex
    peaks.info <- iso.peaks@Info
    nscan <- floor(param.extract@nscan.quant/2)

    smp.classes <- unique(smp.groups[, 'class'])
    names(smp.classes) <- smp.classes
    smpnm.classes <- lapply(smp.classes, function(smp.class) {
      make.names(smp.groups$samplenames[smp.groups[, 'class'] == smp.class])
    })
    nsmp.classes <- sapply(smp.classes, function(smp.class) {
      sum(smp.groups[, 'class'] == smp.class)
    })

    exp.param <- iso.peaks@ExperimentParam
    d.res <- exp.param@d.res
    nms.ft <- names(iso.targets@Isotopologue)
    names(nms.ft) <- nms.ft

    seq.quant <- seq(param.extract@nscan.quant)
    col.nms <- c(seq.quant, 'baseline')
    info.ft <- lapply(nms.ft, function(nm.ft) {
      profile.ft <- peaks.profile[[nm.ft]]
      scans.quant <- sapply(peaks.apex[nm.ft, ], `+`, seq(-nscan, nscan))
      info.mz <- lapply(names(profile.ft), function(nm.mz) {
        profile.mz <- profile.ft[[nm.mz]]
        t(sapply(names(profile.mz), function(nm.smp) {
          profile.smp <- profile.mz[[nm.smp]]
          scans.smp <- scans.quant[, make.names(nm.smp)]

          if (class(profile.smp) == 'data.frame') {
            nscan <- nrow(profile.smp)
            if (any(scans.smp <= 0)) {
              scans.smp <- scans.smp + sum(scans.smp <= 0)
            }
            if (scans.smp >= nscan) {
              scans.smp <- scans.smp - sum(scans.smp >= nscan)
            }
            int <- profile.smp[scans.smp, 'intensity']
            baseline <- peaks.info[[nm.ft]][[nm.smp]][nm.mz, 'baseline']
            res <- c(int, baseline)
          } else {
            res <- rep(NA, length(scans.smp) + 1)
          }
          names(res) <- col.nms
          res
        }))
      })
      names(info.mz) <- names(profile.ft)
      info.mz
    })

    quant.raw <- lapply(nms.ft, function(nm.ft) {
      nms.mz <- names(info.ft[[nm.ft]])
      quant.mz <- t(sapply(nms.mz, function(nm.mz) {
        info.mz <- info.ft[[nm.ft]][[nm.mz]]
        tbl.raw.quant <- info.mz
        tbl.mod.quant <- switch(
          param.extract@method.quant,
          'rawinroi' = {
            tbl.raw.quant[is.na(tbl.raw.quant[, 'baseline']), ] <- NA
            tbl.raw.quant
          },
          'raw' = {
            tbl.raw.quant
          },
          'baseline' = {
            tbl.raw.quant[is.na(tbl.raw.quant[, 'baseline']), ] <- NA
            tbl <- tbl.raw.quant[, seq.quant] - tbl.raw.quant[, 'baseline']
            tbl[tbl < 0] <- 0
            cbind(tbl, 'baseline' = tbl.raw.quant[, 'baseline'])
          })
        rowSums(tbl.raw.quant[, seq.quant])
      }))
      quant.mz <- data.frame('name' = nm.ft,
                             'mz' = as.numeric(rownames(quant.mz)),
                             quant.mz, stringsAsFactors = FALSE)
    })

    quant.table.raw <- do.call(rbind, quant.raw)
    write.csv(quant.table.raw,
              file.path(d.res, 'QuantTableRaw.csv'))

    quant.grps <- lapply(smp.classes, function(smp.class) {
      smpnm.class <- smpnm.classes[[smp.class]]
      quant.grp <- quant.table.raw[, smpnm.class]
      frac.present <- rowSums(!is.na(quant.grp))/nsmp.classes[smp.class]
      quant.grp[frac.present < param.extract@minfrac,
                seq(nsmp.classes[smp.class])] <- NA
      quant.grp
    })

    quant.table.grps <- cbind(quant.table.raw[, 1:2],
                              do.call(cbind, unname(quant.grps)))
    rownames(quant.table.grps) <- NULL
    write.csv(quant.table.grps,
              file.path(d.res, 'QuantTableGroupFilter.csv'))

  })

