# This function maps a matrix of samples (one column per sample) onto a pre-defined "sample-type" space
# Input arguments:
#      samples: The matrix of samples (one column per sample) that we want to map onto the sample-type space,
#   svd_object: The svd of the space to which we are mapping,
#         rank: The rank of the space to which we are mapping (typically, we map onto a relatively low-rank space).

mapToSpace <- function(samples, svd_object, rank) {
  # Allocate memory for the sample space
  sample_space <- matrix(NA, ncol = ncol(samples), nrow = rank)
  # Map each column onto the sample space
  for (i in 1:ncol(samples)) {
    sample_space[, i] <- samples[, i] %*% svd_object$u[, 1:rank]
  }
  # Name the mapped columns in accordance with the original samples
  colnames(sample_space) <- colnames(samples)
  # Return the sample_space matrix
  return(sample_space)
}
