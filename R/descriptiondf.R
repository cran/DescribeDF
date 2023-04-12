#'@title Descriptive Statistics of a Data Frame
#'
#' @param df Data Frame with first column as serial number or date
#' @import dplyr, psych, e1071, stats
#'
#' @return
#' \itemize{
#'   \item desc_df - A table contains 10 descriptive statistics rowwise
#' }
#' @export
#'
#' @examples
#' # create a vector of dates
#'dates <- seq(as.Date("2021-01-01"), as.Date("2021-01-05"), by = "day")

#'# create vectors of random numeric data for columns A through E
#'A <- runif(5, 0, 1)
#'B <- runif(5, 0, 1)
#'C <- runif(5, 0, 1)
#'D <- runif(5, 0, 1)
#'E <- runif(5, 0, 1)

#'# combine the vectors into a data frame
#'df <- data.frame(Date = dates, A = A, B = B, C = C, D = D, E = E)

#'# print the data frame
#'print(df)

#'# Data Description
#'df_descstat(df)
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }

# Descriptive Statistics ####
df_descstat <- function(df) {

  df <- df[, -1, drop = FALSE]  # remove first column

  stats <- c('N', 'Minimum', 'Maximum', 'Mean', 'SD', 'Cond_SD', 'CV(%)',
             'Skewness', 'Kurtosis', 'Shapiro-Wilk')

  res <- matrix(nrow = length(stats), ncol = ncol(df))
  colnames(res) <- colnames(df)
  rownames(res) <- stats

  for (i in 1:ncol(df)) {
    Y <- df[, i]
    sw_test <- shapiro.test(Y)
    p_value <- sw_test$p.value
    sw_stat <- round(sw_test$statistic, 3)

    if (p_value <= 0.01) {
      shapiro_wilk <- "***"
    } else if (p_value <= 0.05) {
      shapiro_wilk <- "**"
    } else if (p_value <= 0.1) {
      shapiro_wilk <- "*"
    } else {
      shapiro_wilk <- ""
    }

    # embedding
    diff.1 <- embed(Y, 2)
    Yt <- diff.1[,1]
    Yt_1 <- diff.1[,2]
    y <- log(Yt/Yt_1)

    res[, i] <- c(length(Y),
                  round(min(Y), 3),
                  round(max(Y), 3),
                  round(mean(Y), 3),
                  round(sd(Y), 3),
                  round(sd(y), 3),
                  round((sd(Y) / mean(Y) * 100), 3),
                  round((psych::skew(Y)), 3),
                  round(e1071::kurtosis(Y), 3),
                  paste0(sw_stat, shapiro_wilk))
  }

  desc_table <- data.frame(Statistics = stats)
  for (i in 1:ncol(df)) {
    desc_table[, colnames(df)[i]] <- res[, i]
  }
  desc_df <- data.frame(desc_table)
  # desc_df_t <-data.frame(t(desc_table))
  #
  # stat_list <- list(desc_df, desc_df_t)
  # return(stat_list)
  return(desc_df)
}

#'@title Non linearity test of a Data Frame
#' @description This function (df_nonlinearity) will give non linearity test result for a df excluding the first column (contains serial number or date). This will give a list of data frames. Data frames are named as the names of columns of the data frame. First column mentions different statistics (eps). Other columns are the Statistics values of the particular dimension.
#' @param df Data Frame with first column as serial number or date
#' @import tseries fNonlinear
#'
#' @return
#' \itemize{
#'   \item result_list - List of data frames named as the column names of provided data frame. Each df is such that first column mentions different statistics and other columns are the Statistics values of the particular dimension.
#' }
#' @export
#'
#' @examples
#' # Create a sequence of numbers from 1 to 100
#'serial <- 1:100

#'# Create six vectors of random numbers, one for each column
#'col1 <- rnorm(100)
#'col2 <- rnorm(100)
#'col3 <- rnorm(100)
#'col4 <- rnorm(100)
#'col5 <- rnorm(100)
#'col6 <- rnorm(100)

# Combine the vectors into a data frame
#'df <- data.frame(serial, col1, col2, col3, col4, col5, col6)
#'df_nonlinearity(df)
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., ... & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }
# Non_Linearity Checking ####

df_nonlinearity <- function(df) {
  result_list <- list()

  for (i in 2:ncol(df)) {
    result <- bdsTest(df[, i])
    statistic <- result@test$statistic
    stats <- data.frame(statistic)
    stats <- round(stats, 3) # round the statistics to 3 decimal places
    p <- data.frame(result@test$p.value)

    for (j in seq_len(nrow(p))) {
      for (k in seq_len(ncol(p))) {
        if (is.na(p[j, k])) {
          stats[j, k] <- paste0(stats[j,], "na")
        } else if (p[j, k] <= 0.01) {
          stats[j, k] <- paste0(stats[j,], "***")
        } else if (p[j, k] <= 0.05) {
          stats[j, k] <- paste0(stats[j,], "**")
        } else if (p[j, k] <= 0.1) {
          stats[j, k] <- paste0(stats[j,], "*")
        } else {
          stats[j, k] <- stats[j,]
        }
      }
    }

    result_df <- data.frame(cbind(c('eps[1]', 'eps[2]', 'eps[3]', 'eps[4]'),
                                  c(stats[1,], stats[3,], stats[5,], stats[7,]),
                                  c(stats[2,], stats[4,], stats[6,], stats[8,])))

    colnames(result_df) <- c('Statistics', 'm=2', 'm=3')

    col_name <- colnames(df)[i]

    result_list[[col_name]] <- result_df
  }

  return(result_list)
}

#'@title Stationarity tests of a Data Frame
#' @description This function (df_stationarity) will give a list of three data frames: 'ADF', 'PP', 'KPSS'. This will also indicate whether the data is stationary or not according to the null hypothesis of the corresponding tests. The data frame must contain serial number or date or anything in the 1st column.This function will exclude the 1st column of the data frame and will perform tests on other columns.
#' @param df Data Frame with first column as serial number or date
#' @import tseries
#'
#' @return
#' \itemize{
#'   \item test_results - List of three data frames: 'ADF', 'PP', 'KPSS'
#' }
#' @export
#'
#' @examples
#' # create a vector of dates
#'dates <- seq(as.Date("2021-01-01"), as.Date("2021-01-05"), by = "day")

#'# create vectors of random numeric data for columns A through E
#'A <- runif(5, 0, 1)
#'B <- runif(5, 0, 1)
#'C <- runif(5, 0, 1)
#'D <- runif(5, 0, 1)
#'E <- runif(5, 0, 1)

#'# combine the vectors into a data frame
#'df <- data.frame(Date = dates, A = A, B = B, C = C, D = D, E = E)

#'# print the data frame
#'print(df)

#'# stationarity results
#'df_stationarity(df)
#' @references
#' \itemize{
#'\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.

#' \item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., ... & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.

#' }

# Stationarity Checking ####

df_stationarity <- function(df) {
  # This function will give a list of three data frames: 'ADF', 'PP', 'KPSS'
  # This function will also indicate whether the data is stationary or not.
  # The data frame must contain serieal number or date or anything in the 1st column
  # This function will exclude the 1st column of the data frame and will perform tests on other columns
  # Get column names
  cols <- names(df)[-1]

  # Create empty data frames for results
  adf_df <- data.frame(Test = "ADF",
                       Parameters = c("Dickey-Fuller", "Lag order", "p-value", "Stationary?"),
                       stringsAsFactors = FALSE)
  pp_df <- data.frame(Test = "PP",
                      Parameters = c("Dickey-Fuller (alpha)", "Truncation Lag parameter", "p-value", "Stationary?"),
                      stringsAsFactors = FALSE)
  kpss_df <- data.frame(Test = "KPSS",
                        Parameters = c("KPSS level", "Truncation Lag parameter", "p-value", "Stationary?"),
                        stringsAsFactors = FALSE)

  # Apply tests to each column and store results in data frames
  for (col in cols) {
    adf_res <- suppressWarnings(adf.test(df[, col]))
    adf_stationary <- ifelse(adf_res$p.value <= 0.01, "Yes***", "No***")
    adf_df[[col]] <- c(round(adf_res$statistic, 3), adf_res$parameter, round(adf_res$p.value, 3), adf_stationary)

    pp_res <- suppressWarnings(pp.test(df[, col]))
    pp_stationary <- ifelse(pp_res$p.value <= 0.01, "Yes***", "No***")
    pp_df[[col]] <- c(round(pp_res$statistic, 3), pp_res$parameter, round(pp_res$p.value, 3), pp_stationary)

    kpss_res <- suppressWarnings(kpss.test(df[, col]))
    kpss_stationary <- ifelse(kpss_res$p.value >= 0.01, "Yes***", "No***")
    kpss_df[[col]] <- c(round(kpss_res$statistic, 3), kpss_res$parameter, round(kpss_res$p.value, 3), kpss_stationary)
  }

  # Combine data frames into a list
  test_results <- list(ADF = data.frame(adf_df), PP = data.frame(pp_df), KPSS = data.frame(kpss_df))

  return(test_results)
}

