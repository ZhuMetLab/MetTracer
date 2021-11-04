##############################
## Parameter classes
##############################
## BasicParam class
setClass("ParamMetTracer",
         representation = representation("VIRTUAL"))

setClassUnion("nullOrCharacter", c("NULL", "character"))

## Experiment setup
setClass("ExperimentParam",
         slots = c(
           wd = "character",
           d.res = 'character',
           d.tmp = 'character',
           element.trace = 'character',
           element.label = 'character',
           equipment = 'character',
           resolution = 'numeric',
           ppm = 'numeric',
           res.define = 'numeric',
           nSlaves = 'numeric'
         ))

## generating extractTargets parameters
setClass('IsotopologueParam',
         slots = c(file.id = 'character',
                   file.xset = 'character',
                   method = 'character',
                   value = 'character',
                   rt.extend = 'numeric',
                   is.plot.pattern = 'logical'),
         contains = c("ParamMetTracer"),
         prototype = prototype(
           file.id = '../unlabelled/MetProcess-Result/MRN_annotation_result/MRN.annotation.result.csv',
           file.xset = '../unlabelled/xset-centWave-final.Rda',
           method = 'max',
           value = 'into',
           rt.extend = 30,
           is.plot.pattern = FALSE),
         validity = function(object) {
           msg <- character()
           if (length(object@file.id) != 1)
             msg <- c(msg, paste0("'file.id' has to be character",
                                  " of length 1."))
           if (length(object@file.xset) != 1)
             msg <- c(msg, paste0("'file.xset' has to be character",
                                  " of length 1."))
           if (length(object@method) != 1)
             msg <- c(msg, paste0("'method' has to be character",
                                  " of length 1."))
           if (length(object@value) != 1)
             msg <- c(msg, paste0("'value' has to be character",
                                  " of length 1."))
           if (length(object@rt.extend) != 1 || object@rt.extend < 0)
             msg <- c(msg, paste0("'rt.extend' has to be numeric >= 0",
                                  " of length 1."))
           if (length(object@is.plot.pattern) != 1)
             msg <- c(msg, paste0("'is.plot.pattern' has to be logical",
                                  " of length 1."))
           if (length(msg))
             msg
           else
             TRUE
         }
)

setClass('ExtractParam',
         slots = c(
           d.extract = 'character',
           range.extract = 'character',
           method.align = 'character',
           method.quant = 'character',
           method.best = 'character',
           nscan.quant = 'numeric',
           minfrac = 'numeric',
           adj.unlabel = 'nullOrCharacter',
           adj.label = 'nullOrCharacter',
           adj.contaminate = "logical",
           thr.contaminate = "numeric"),
         contains = c("ParamMetTracer"),
         prototype = prototype(
           d.extract = 'extract',
           range.extract = 'all',
           method.align = 'apex',
           method.best = 'maxint',
           method.quant = 'rawinroi',
           nscan.quant = 3,
           minfrac = 0.5,
           adj.unlabel = NULL,
           adj.label = NULL,
           adj.contaminate = FALSE,
           thr.contaminate = 0.02),
         validity = function(object) {
           msg <- character()
           if (length(object@d.extract) != 1)
             msg <- c(msg, paste0("'d.extract' has to be character",
                                  " of length 1."))
           if (length(object@range.extract) != 1 || object@range.extract < 0)
             msg <- c(msg, paste0("'range.extract' has to be character",
                                  " of length 1."))
           if (length(object@method.align) != 1)
             msg <- c(msg, paste0("'method.align' has to be character",
                                  " of length 1."))
           if (length(object@method.best) != 1)
             msg <- c(msg, paste0("'method.best' has to be character",
                                  " of length 1."))
           if (length(object@method.quant) != 1)
             msg <- c(msg, paste0("'method.quant' has to be character",
                                  " of length 1."))
           if (length(object@nscan.quant) != 1)
             msg <- c(msg, paste0("'nscan.quant' has to be numeric",
                                  " of length 1."))
           if (length(object@minfrac) != 1 | object@minfrac < 0 | object@minfrac > 1)
             msg <- c(msg, paste0("'minfrac' has to be numeric in [0, 1]",
                                  " of length 1."))
           if (length(object@thr.contaminate) != 1)
             msg <- c(msg, paste0("'thr.contaminate' has to be numeric",
                                  " of length 1."))
           if (length(msg))
             msg
           else
             TRUE
         }
)

setClass('PeakdetectionParam',
         slots = c(peakwidth = 'numeric',
                   snthr = 'numeric',
                   prefilter = 'numeric',
                   fitgauss = 'logical',
                   smoothNoise = 'numeric',
                   method.peakdetection = 'character',
                   method.roi = 'character',
                   method.smooth = 'character',
                   method.baseline = 'character'),
         contains = c("ParamMetTracer"),
         prototype = prototype(
           peakwidth = c(5, 30),
           snthr = 3,
           prefilter = c(3, 100),
           fitgauss = FALSE,
           smoothNoise = 0.35,
           method.peakdetection = 'centWave',
           method.roi = 'continuous',
           method.smooth = 'gaussian',
           method.baseline = 'xcms'
         ),
         validity = function(object) {
           msg <- character()
           if (length(object@peakwidth) != 2)
             msg <- c(msg, paste0("'peakwidth' has to be numeric",
                                  " of length 2."))
           if (length(object@snthr) != 1 || object@smoothNoise < 0)
             msg <- c(msg, paste0("'snthr' has to be numeric >= 0",
                                  " of length 1."))
           if (length(object@prefilter) != 2)
             msg <- c(msg, paste0("'prefilter' has to be numeric",
                                  " of length 2."))
           if (length(object@fitgauss) != 1)
             msg <- c(msg, paste0("'fitgauss' has to be logical",
                                  " of length 1."))
           if (length(object@smoothNoise) != 1 || object@smoothNoise < 0)
             msg <- c(msg, paste0("'smoothNoise' has to be numeric >= 0",
                                  " of length 1."))
           if (length(object@method.peakdetection) != 1)
             msg <- c(msg, paste0("'method.peakdetection' has to be character",
                                  " of length 1."))
           if (length(object@method.roi) != 1)
             msg <- c(msg, paste0("'method.roi' has to be character",
                                  " of length 1."))
           if (length(object@method.smooth) != 1)
             msg <- c(msg, paste0("'method.smooth' has to be character",
                                  " of length 1."))
           if (length(object@method.baseline) != 1)
             msg <- c(msg, paste0("'method.baseline' has to be character",
                                  " of length 1."))
           if (length(msg))
             msg
           else
             TRUE
         }
)

##############################
## Data classes
##############################

setClass('IsotopologueTargets',
         slots = c(TargetInfo = 'data.frame',
                   Isotopologue = 'list',
                   Formula = 'list',
                   Profile = 'list',
                   rtRange = 'matrix',
                   mzRange = 'list',
                   mzIsomer = 'list',
                   ExperimentParam = 'ExperimentParam',
                   .processHistory = 'list'))

setGeneric('GenerateIsotopologue',  function(param, experimentParam, ...)
  standardGeneric('GenerateIsotopologue'))

setClass('IsotopicPeaks',
         slots = c(TargetInfo = 'data.frame',
                   SampleGroups = 'data.frame',
                   Peaks = 'data.frame',
                   Profile = 'list',
                   Apex = 'data.frame',
                   Info = 'list',
                   ExperimentParam = 'ExperimentParam',
                   .processHistory = 'list'))

setClass('IsotopologueExtracted',
         slots = c(TargetInfo = 'data.frame',
                   Isotopologue = 'list',
                   Profile = 'list',
                   ExperimentParam = 'ExperimentParam',
                   .processHistory = 'list'))


setGeneric('ExtractIsotopicPeaks',  function(param.extract, param.pd, iso.targets, ...)
  standardGeneric('ExtractIsotopicPeaks'))

setGeneric('ExtractIsotopologue',  function(param.extract, iso.peaks, iso.targets, ...)
  standardGeneric('ExtractIsotopologue'))

setMethod('show', 'IsotopologueExtracted', function(object) {
  cat("IsotopologueExtracted object of", nrow(object@TargetInfo), "features:\n")
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/3^10, 3), "MB\n")
})

setMethod('show', 'IsotopicPeaks', function(object) {
  cat("IsotopicPeaks object of", nrow(object@TargetInfo), "features:\n")
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/3^10, 3), "MB\n")
})

setMethod('show', 'IsotopicPeaks', function(object) {
  cat("IsotopicPeaks object of", nrow(object@TargetInfo), "features:\n")
  memsize <- object.size(object)
  cat("Memory usage:", signif(memsize/3^10, 3), "MB\n")
})
