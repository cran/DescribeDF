\name{df_descstat}
\alias{df_descstat}
\title{Descriptive Statistics of A Data Frame}
\usage{
df_descstat(df)
}
\arguments{
\item{df}{Data Frame with first column as serial number or date}
}
\value{
\itemize{
  \item desc_df - A table contains 10 descriptive statistics rowwise
}
}
\description{
Gives descriptive statistics for a data frame (df) excluding the first column. First column mentions different statistics. Other columns are the Statistics values of the particular column (series). "p_value <= 0.01: ***; p_value <= 0.05: **; p_value <= 0.1: *".
}
\examples{
# create a vector of dates
dates <- seq(as.Date("2021-01-01"), as.Date("2021-01-05"), by = "day")

# create vectors of random numeric data for columns A through E
A <- runif(5, 0, 1)
B <- runif(5, 0, 1)
C <- runif(5, 0, 1)
D <- runif(5, 0, 1)
E <- runif(5, 0, 1)

# combine the vectors into a data frame
df <- data.frame(Date = dates, A = A, B = B, C = C, D = D, E = E)

# print the data frame
print(df)

# Data Description
df_descstat(df)
}
\references{
\itemize{
\item Garai, S., & Paul, R. K. (2023). Development of MCS based-ensemble models using CEEMDAN decomposition and machine intelligence. Intelligent Systems with Applications, 18, 200202.
\item Garai, S., Paul, R. K., Rakshit, D., Yeasin, M., Paul, A. K., Roy, H. S., Barman, S. & Manjunatha, B. (2023). An MRA Based MLR Model for Forecasting Indian Annual Rainfall Using Large Scale Climate Indices. International Journal of Environment and Climate Change, 13(5), 137-150.
}
}

