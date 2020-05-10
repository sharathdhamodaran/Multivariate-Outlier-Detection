boxplot_outliers <- function(dt,
                             input_b = 4,
                             iqr_mult = 1.5) {

  # Compute adjusted rule based on ‘medcouple’, a robust estimator of skewness

  # Args:
  #   dt: datatable
  #   input_b: a constant
  #   iqr_mult: Inter Quartile Range Multiplier
  #
  # Returns:
  #   A List of rules applied to the dataframe

  my_mc <- robustbase::mc(dt)
  qa <- quantile(dt, probs = c(.25, .75))
  iqr <- qa[2] - qa[1]
  standard_rule <- qa[2] + iqr_mult * iqr
  adjusted_rule <- qa[2] + iqr_mult * exp(input_b * my_mc) * iqr
  return(list(standard_rule = standard_rule, adjusted_rule = adjusted_rule))
}
