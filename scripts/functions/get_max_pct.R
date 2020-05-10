get_max_pct <- function(x) {

  # Get the maximum percentage of na's accross columns

  # Args:
  #   x: vector
  #
  # Returns:
  #   Max percentage of na's for the vector

  na_count <- sum(is.na(x))

  d <- data.frame(value = x)
  data.table::setDT(d)
  val <- d[, .(count = .N), by = value]

  if (nrow(val) == 0) {
    return(1)
  } else {
    if (na_count > 0) {
      return(max(na_count, val$count) / length(x))
    } else {
      return(max(val$count) / length(x))
    }
  }
}
