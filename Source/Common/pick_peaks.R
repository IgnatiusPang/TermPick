# Spike detection above background ----------------------------------------
# Reference for the codes: https://github.com/cran/peakPick/blob/master/R/peakpicking.R

#' Helper functions for computing scaled peaks
#'
#' @param mat matrix of series with series organized columnwise
#' @param positions vector of positions to exclude from threshold computation.
#'   Internal use only; follows R's rules of matrix indexing by vector.
#' @param nsd numeric number of standard deviations for limits (see Value)
#' @return vector limits defined as means + nsd SEMs computed for the columns of
#'   mat, excluding positions from the calculation.
#' @keywords internal
#'
scaledrow <- function (mat, positions, nsd) {
  # assume matrix with vector of positions to exclude when computing sems.
  # calculate values scaled by sems for center row.
  # center determined by winlen + 1.
  matexclude <- mat
  matexclude[positions] <- NA
  means <- colMeans(matexclude, na.rm=TRUE)
  sds <- matrixStats::colSds(matexclude, na.rm=TRUE)
  # scale by sqrt(n)
  sems <- sds/sqrt(matrixStats::colCounts(!is.na(matexclude)))
  return(means + nsd*sems)
}

#' Detects spikes in series of numbers.
#'
#' This algorithm detects spikes rising above a user-specified number of
#' standard deviations numbers in a certain window.  Use this algorithm to
#' detect short spikes rather than smooth bumps in series of numbers.  Please
#' refer to the paper by Weber et al. for more details.
#'
#' @param mat matrix of series with series organized columnwise.  The algorithm
#'   treats each column separately.
#' @param roi vector of two integers (c(min, max)) defining positions in all
#'   series (rows in mat) to consider for spike detection, used together with
#'   winlen.  Must lie within the interval [2, nrow(mat) - 1].  Will be coerced
#'   to integers.
#' @param winlen integer defining the window of positions to consider for mean
#'   and sem estimation for each series.  Each estimation limits itself to the
#'   position and a plus/minus winlen positions large window.  Thus, winlen must
#'   not be chosen larger than that the windows fit within mat, given the roi.
#'   I.e. roi[1] - winlen >=1 AND roi[length(roi)] + winlen <= nrow(mat).  Will
#'   be coerced to an integer.
#' @param spike.min.sd numeric minimum number of standard deviations for a spike
#'   to rise above the mean in order to be considered for a spike call and to be
#'   excluded from the mean estimation of each subsequent iteration of the spike
#'   calling algorithm
#' @param mc.cores the number of cores do perform this calculation
#' @param verbose Boolean indicating the number of new peaks detected with each
#'   iteration.  The algorithm stops as soon as this number does not sink
#'   anymore.  Turn this on if running into problems.
#' @return boolean matrix corresponding to mat, representing spike positions.
#' @references Weber, C.M., Ramachandran, S., and Henikoff, S. (2014).
#'   Nucleosomes are context-specific, H2A.Z-modulated barriers to RNA
#'   polymerase. Molecular Cell 53, 819-830.
#' @export
detect.spikes <- function (mat, roi, winlen, spike.min.sd=3,
                           #mc.cores=1,
                           verbose=FALSE) {

  ## check matrix validity
  if (!is.matrix(mat) || nrow(mat) < 3L)
    stop("mat must be a matrix with at least 3 rows")

  ## check ROI validity
  if (!is.integer(winlen))
    winlen <- as.integer(winlen[1])
  if (!is.integer(roi))
    roi <- as.integer(roi)
  if (length(roi) != 2L || diff(roi) < 1L)
    stop("roi must contain exatly two elements, c(min, max), with max > min")
  if (roi[1] < 2L || roi[2] > nrow(mat) - 1L)
    stop("roi must lie within the interval [2, nrow(mat) - 1]")
  if (roi[1] - winlen < 1L || roi[2] + winlen > nrow(mat))
    stop("winlen too large for the roi")

  ## extract ROI
  roivec <- roi[1]:roi[2]
  roimat <- mat[roivec, ]
  ## build logical matrix, same dimensions as mat, will record spikes
  flaggedstall <- matrix(FALSE, nrow(mat), ncol(mat))
  colnames(flaggedstall) <- colnames(mat)

  ## iterative algorithm:
  changed <- integer(0)
  repeat {
    ## successively compute for all rows in roi a threshold mean + 3*sem.
    ## Exclude already detected peaks from mean and sem estimation (second
    ## argument to scaledrow)
    test <- do.call(rbind, furrr::future_map(roivec, function (rownum) {
      scaledrow(mat[(rownum-winlen):(rownum+winlen), ],
                which(flaggedstall[(rownum-winlen):(rownum+winlen), ]),
                nsd=spike.min.sd)
    })) # , mc.cores=mc.cores

    ## new iteration of spike hits
    newflaggedstall <- flaggedstall
    newflaggedstall[roivec, ] <- roimat > test

    ## if new iteration same as previous, stop
    if (identical(newflaggedstall, flaggedstall))
      return(flaggedstall)

    ## record how many new peaks detected
    changed <- append(changed, length(which(newflaggedstall != flaggedstall)))
    if (verbose)
      print(changed[length(changed)])

    ## if number of new peaks equal or more than last iteration, we reached a
    ## local minimum; stop
    if (length(changed) > 1 && diff(changed)[length(changed) - 1] >= 0)
      return(flaggedstall)

    ## otherwise, store result of iteration and start over
    flaggedstall <- newflaggedstall
  }
}
