get_params_to_remove <- function(df,
                                 max_missing_pct = 0.15,
                                 min_number_of_levels = 2,
                                 max_pct_one_value = 0.5) {

  # Get the numeric columns to remove before Phase I of outlier detection

  # Args:
  #   df: A data frame or data table for cleaning
  #   max_missing_pct: Maximum allowed % of missing values for every column
  #   min_number_of_levels: Minimum number of levels allowed for every column
  #   max_pct_one_value: Maximum % allowed for a single value for every column
  #
  # Returns:
  #   A character vector containing the names of columns that need to be removed
  #   from the input data frame or data table

  if (!data.table::is.data.table(df)) data.table::setDT(df)

  df <- df[, sapply(df, is.numeric), with = FALSE]

  # 1. Columns with greater than max_missing_pct missing values are removed
  # 2. Columns that have <= min_number_of_levels levels are removed
  # 3. Columns with one value > max_pct_one_value of all rows are removed
  # 4. If the correlation between 2 columns is 1, one of the columns is removed

  test1 <- apply(df, 2, function(x) sum(is.na(x)) / length(x))
  test2 <- apply(df, 2, function(x) length(unique(x)))

  test3 <- apply(df, 2, get_max_pct)
  test <- data.frame(cbind(test1, test2, test3))
  columns_to_remove_1 <- rownames(test[(test$test1 > max_missing_pct |
    test$test2 <= min_number_of_levels |
    test$test3 > max_pct_one_value), ])

  df[, (columns_to_remove_1) := NULL]
  cor.mat <- cor(df)
  cor.1 <- data.frame(which(cor.mat >= 0.99999999, arr.ind = TRUE))
  r <- cor.1[cor.1$row != cor.1$col, ]

  var1 <- rownames(cor.mat)[r[, 1]]
  var2 <- colnames(cor.mat)[r[, 2]]
  test4.temp <- data.frame(cbind(var1, var2))
  test4 <- head(test4.temp, nrow(test4.temp) / 2)

  columns_to_remove_4 <- as.character(test4$var1)

  columns_to_remove <- unique(c(columns_to_remove_1, columns_to_remove_4))

  return(columns_to_remove)
}
