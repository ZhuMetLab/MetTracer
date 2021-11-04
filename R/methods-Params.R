#' Experiment setup
#'
#' @param wd Working directory for experiment to be processed
#' @param d.tmp Folder for storing files when processing. By default, it is the
#'      same with \code{wd}
#' @param d.res Folder for storing results
#' @param element.trace Monoisotope elements, e.g., 12C or 14N
#' @param element.label Isotope labelled elements, e.g., 13C or 15N
#' @param equipment Equipment type: QTOF/Obiwarp
#' @param resolution Resolution of equipment, for QTOF, it is 3000 by default
#' @param ppm Mass to charge value (m/z) tolerance for MS1 ions
#' @param res.define Where the PPP resolution starts from. For QTOF equipments
#'      in zhulab, it is set to 400
#' @param nSlaves Number of threads to be used
#' @return an \code{ExperimentParam} object
#' @export
ExperimentParam <- function(
  wd = ".",
  d.res = file.path(wd, 'Result'),
  d.tmp = file.path(d.res, 'tmp'),
  element.trace = '12C',
  element.label = '13C',
  equipment = c('QTOF', 'Obiwarp'),
  resolution = NULL,
  ppm = 25,
  res.define = 400,
  nSlaves = 4) {
  options(mc.cores = nSlaves)

  if (!dir.exists(wd)) {
    stop(paste("Working directory setup wrong!\n The file path",
               wd, "is not found!"))
  }
  if (!dir.exists(d.res)) {
    dir.create(d.res, recursive = TRUE)
  }
  if (!dir.exists(d.tmp)) {
    dir.create(d.tmp, recursive = TRUE)
  }

  equipment <- match.arg(equipment)
  if (missing(resolution) | is.null(resolution)) {
    resolution <- switch(equipment, 'QTOF' = 3000, 'Obiwarp' = 10000)
  }

  return(new("ExperimentParam",
             wd = wd,
             d.res = d.res,
             d.tmp = d.tmp,
             element.trace = element.trace,
             element.label = element.label,
             equipment = equipment,
             resolution = resolution,
             ppm = ppm,
             res.define = res.define,
             nSlaves = nSlaves))
}


#' Isotopoluage target generation parameter setup
#'
#' @param file.id MetDNA annotation result file
#' @param file.xset xcms peak detection file
#' @param method method for creating representative peak profile, (max, sum, median)
#' @param value value field for creating representative peak profile (into, intb, maxo)
#' @param rt.extend extention thresholds for RT
#' @param is.plot.pattern if plot isotopoluage patterns
#' @return an \code{extractTargetsParam} object
#' @export
IsotopologueParam <- function(
  file.id = 'unlabelled/MetProcess-Result/MRN_annotation_result/MRN.annotation.result.csv',
  file.xset = 'unlabelled/xset-centWave-final.Rda',
  method = c('max', 'sum', 'median'),
  value = c('into', 'intb', 'maxo'),
  rt.extend = 30,
  is.plot.pattern = FALSE) {

  method <- match.arg(method)
  value <- match.arg(value)

  return(new("IsotopologueParam",
             file.id = file.id,
             file.xset = file.xset,
             method = method,
             value = value,
             rt.extend = rt.extend,
             is.plot.pattern = is.plot.pattern))
}


#' Isotopoluage extraction parameter setup
#'
#' @param d.extract directory of samples to be extracted
#' @param range.extract EIC extraction range
#' \itemize{
#'     \item[] 'all' - extract the entire eic
#'     \item[] 'roi' - extract the rois related to rt range of target peaks
#'     \item []'rtrange' - extract the eic at the range of rt range of target peaks
#'     }
#' @param method.align method for creating representative peak profile
#' \itemize{
#'     \item[] 'apex' align with the apex of the most similar peaks
#'     \item[] 'eic' align with the correlation of the eics from unlabelled and
#'           labelled samples
#' }
#' @param method.best method for finding best correlated peak groups for extracted
#'     isotopoluage to determine the rt range of extracted peaks/eics
#'  \itemize{
#'      \item[] 'maxint' select the highest peak as rt range reference
#'      \item[] 'maxcor select the most similar peak to representive peak as rt range reference
#'  }
#'  @param method.quant method for quantifying extracted isotopic peaks.
#'  \itemize{
#'      \item[] 'rawinroi' sum with raw intensities of significant scans (in ROI)
#'      \item[] 'baseline' baseline substraxted rawinroi method
#'      \item[] 'raw' sum with raw intensities anyway
#'  }
#' @param nscan.quant number of scans to sum up for quantifying
#' @param minfrac minimal fraction threshold to determine the isotopic peak
#'      existed or not in the sample groups
#' @param adj.unlabel character vector, group names of unlabelled samples for
#'  adjusting the isotopoluage of labelled samples
#' @param adj.label character vector, group names of labelled samples to be
#'  adjusted using the unlabelled samples (pairwised with the sample groups in
#'  \code{adj.unlabel})
#' @param adj.contaminate logical: if adjust the contaminate using unlabelled
#'  samples, if \code{TRUE}, \code{adj.unlabel} and \code{adj.label} must be set
#'
#' @return an \code{ExtractParam} object
#' @export
ExtractParam <- function(
  d.extract = 'extract',
  range.extract = c('all', 'roi', 'rtrange'),
  method.align = c('apex', 'eic'),
  method.best = c('maxint', 'maxcor'),
  method.quant = c('rawinroi', 'raw', 'baseline'),
  nscan.quant = 3,
  minfrac = 0.5,
  adj.unlabel = NULL,
  adj.label = NULL,
  adj.contaminate = FALSE,
  thr.contaminate = 0.02) {

  range.extract <- match.arg(range.extract)
  method.align <- match.arg(method.align)
  method.best <- match.arg(method.best)
  method.quant <- match.arg(method.quant)

  if (adj.contaminate) {
    if (is.null(adj.label)) {
      stop("The following parameters must be set first: \n adj.unlagel and adj.label")
    }
  }

  if (length(adj.unlabel) != length(adj.label)) {
    stop("Length of the following two parameters must be equal:\n adj.unlagel and adj.label")
  }

  return(new("ExtractParam",
             d.extract = d.extract,
             range.extract = range.extract,
             method.align = method.align,
             method.best = method.best,
             method.quant = method.quant,
             nscan.quant = nscan.quant,
             minfrac = minfrac,
             adj.unlabel = adj.unlabel,
             adj.label = adj.label,
             adj.contaminate = adj.contaminate,
             thr.contaminate = thr.contaminate))
}

#' Peak detection parameter setup
#'
#' @param peakwidth `numeric(2)` with the lower and upper bound of the
#'     expected peak width.
#' @param snthr `numeric(1)` defining the signal to noise ratio cutoff.
#'     Peaks with a signal to noise ratio < `snthr` are omitted.
#' @param prefilter `numeric(2)` (`c(k, I)`): only regions of interest with at
#'     least `k` centroids with signal `>= I` are returned in the first
#'     step.#'
#' @param fitgauss `logical(1)` whether or not a Gaussian should be fitted
#'     to each peak.
#' @param method.peakdetection `character(1)` Peak detection methods
#' \itemize{
#'     \item[] 'centWave' - find peaks with 'centWave' algorim
#'     \item[] 'localmax' - find peaks by detecting the local max and local min
#' }
#' @param method.roi method for finding rois in extracted EICs.
#' \itemize{
#'     \item[] 'continuous' - traditional centwave roi finding method
#'     \item[] 'aroundMax' - finding local maximums and determing a scan range
#'      based on the peakwidth to define the roi range
#' }
#' @param method.baseline method for determing baselines, only xcms supported currently
#' @param method.smooth `character(1)` method for smoothing the EICs (Gaussian, LOESS and SG)
#' @param method.peakdetection `character(1)` method for peak detection

#' @return an \code{PeakdetectionParam} object
#' @export
PeakdetectionParam <- function(
  peakwidth = c(5, 30),
  snthr = 3,
  prefilter = c(3, 100),
  fitgauss = FALSE,
  method.peakdetection = c('centWave', 'localMax'),
  method.roi = c('continuous', 'aroundMax'),
  method.baseline = c('xcms', 'centWaveP'),
  method.smooth = c('Gaussian', 'LOESS', 'SG')) {

  method.peakdetection <- match.arg(method.peakdetection)
  method.roi <- match.arg(method.roi)
  method.baseline <- match.arg(method.baseline)
  method.smooth <- match.arg(method.smooth)

  return(new("PeakdetectionParam",
             peakwidth = peakwidth,
             snthr = snthr,
             prefilter = prefilter,
             fitgauss = fitgauss,
             method.peakdetection = method.peakdetection,
             method.roi = method.roi,
             method.baseline = method.baseline,
             method.smooth = method.smooth))
}
