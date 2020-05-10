unknown_dec_rule <- function(batch_size,
                             t2_threshold_quantile,
                             t2_outlier_quantile,
                             alpha,
                             beta) {

  # A function to decide if a group is Out Of Control

  # Args:
  #    batch_size: Sample size of the group
  #    t2_threshold_quantile: The alpha value used in Phase I for UCL
  #    t2_outlier_quantile: The % of outlier removed from Phase I UCL
  #    alpha: Type I error
  #    beta: Type II error

  # Returns:
  #    dec_rule: A list containing a decision rule
  #    approximate_alpha: actual type I error
  #    approximate_p1_power: actual type II error

  # Formula Values:
  p0 <- t2_threshold_quantile + t2_outlier_quantile - (t2_threshold_quantile *
    t2_outlier_quantile)
  N <- batch_size
  Za <- qnorm(1 - alpha, 0, 1)
  # Zb <- qnorm(beta, 0, 1)
  # large sample rule:
  dec_rule <- ceiling(N * (Za * sqrt(p0 * (1 - p0) / N) + p0))
  # dec_rule_pct <- 100 * dec_rule / N
  # calculate error rates and effect size:
  index <- 0:dec_rule
  approximate_alpha <- 1 - sum(dbinom(index, size = N, prob = p0))
  effect_size <- pwr::pwr.p.test(
    sig.level = approximate_alpha,
    power = 1 - beta,
    n = N,
    alternative = "greater"
  )$h
  approximate_p1 <- sin((effect_size + 2 * asin(sqrt(p0))) / 2)^2
  approximate_p1_power <- 1 - sum(dbinom(index, size = N, prob = approximate_p1))

  return(list(dec_rule, approximate_alpha, approximate_p1_power))
}
