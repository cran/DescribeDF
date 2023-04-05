\name{df_nonlinearity}
\alias{df_nonlinearity}
\title{Non linearity test of A Data Frame}
\usage{
df_nonlinearity(df)
}
\arguments{
\item{df}{Data Frame with first column as serial number or date}
}
\value{
\itemize{
  \item result_list - List of data frames named as the column names of provided data frame. Each df is such that first column mentions different statistics and other columns are the Statistics values of the particular dimension.
}
}
\description{
Gives non linearity test result for a df excluding the first column (contains serial number or date). This will give a list of data frames. Data frames are named as the names of columns of the data frame. First column mentions different statistics (eps). Other columns are the Statistics values of the particular dimension. "p_value <= 0.01: ***; p_value <= 0.05: **; p_value <= 0.1: *".
}
\examples{
# Create a sequence of numbers from 1 to 100
serial <- 1:100

# Create six vectors of random numbers, one for each column
col1 <- rnorm(100)
col2 <- rnorm(100)
col3 <- rnorm(100)
col4 <- rnorm(100)
col5 <- rnorm(100)
col6 <- rnorm(100)

# Combine the vectors into a data frame
df <- data.frame(serial, col1, col2, col3, col4, col5, col6)
df_nonlinearity(df)
}
\references{
\itemize{
\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.
\item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
}
}

