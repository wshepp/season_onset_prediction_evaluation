library(ggplot2)
library(mgcv)

source("preprocessing.R")

generate_model <- function(cleaned_df){
  # read in data for current plot
  data_df <- data.frame(
    week = cleaned_df$week,
    ilip = cleaned_df$ilip
  )
  max_weeks <- length(cleaned_df$week)
  # generate model for each indicator
  model <- gam(log(ilip+EPS)~s(as.numeric(week), k = as.integer(floor(max_weeks/5))), 
               family = gaussian(), 
               data = data_df)
  
  nsim <- 1000
  # fit simulations for ili+
  pdat <- with(data_df, 
               data.frame(week = week, ilip = log(ilip)))
  sim <- mvrnorm(nsim, mu = coef(model), Sigma = vcov(model))
  lp <- predict(model, newdata = pdat, type = "lpmatrix")
  fits <- lp %*% t(sim)
  
  # find quantiles for each indicator
  plotting_df <- data.frame()
  for (j in 1:max_weeks){
    quan<- quantile(fits[j,], c(0.5,0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975))
    row <- data.frame(time = pdat$week[j],
                      t_step = j,
                      y=exp(quan[[1]]),
                      lb_50 = exp(quan[[5]]),
                      ub_50 = exp(quan[[6]]),
                      lb_80 = exp(quan[[4]]),
                      ub_80 = exp(quan[[7]]),
                      lb_90 = exp(quan[[3]]),
                      ub_90 = exp(quan[[8]]),
                      lb_95 = exp(quan[[2]]),
                      ub_95 = exp(quan[[9]]))
    plotting_df <- rbind(plotting_df, row)
  }
  
  # find growth rate for ili+ sims
  pdat <- with(data_df, 
               data.frame(week = week, ilip = log(ilip)))
  X0 <- predict(model, pdat, type = "lpmatrix")
  pdat <- pdat + 7
  X1 <- predict(model, pdat, type = "lpmatrix")
  Xp <- (X1 - X0)
  df_fits <- Xp %*% t(sim)
  
  plotting_gr_df <- data.frame()
  # find quantiles for each indicator growth rate
  for (j in 1:max_weeks){
    quan<- quantile(df_fits[j,], c(0.5,0.025,0.05,0.1,0.25,0.75,0.9,0.95,0.975))
    row <- data.frame(time = pdat$week[j],
                      t_step = j,
                      y=(quan[[1]]),
                      lb_50 = (quan[[5]]),
                      ub_50 = (quan[[6]]),
                      lb_80 = (quan[[4]]),
                      ub_80 = (quan[[7]]),
                      lb_90 = (quan[[3]]),
                      ub_90 = (quan[[8]]),
                      lb_95 = (quan[[2]]),
                      ub_95 = (quan[[9]]))
    plotting_gr_df <- rbind(plotting_gr_df, row)
  } 
  # plot current graph
  
  return(plotting_gr_df)
}
  
# read in data
flunet_df <- read.csv('Data/VIW_FNT.csv', header = T)
fluid_df <- read.csv('Data/VIW_FID_EPI.csv', header = T)

start_date <- '2012-01-01'
end_date <- '2024-06-30'
country <- 'Australia'

cleaned_df <- preprocessing(country, fluid_df, flunet_df, start_date, end_date)

gr_results <- generate_model(cleaned_df)
gr_results

ggplot(data = gr_results) +
  geom_line(aes(x = time, y = y, color = "Fitted Model")) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = "Fitted Model"
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = "Fitted Model"
  ), alpha = 0.2) +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-07-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = c("Fitted Model" = "black")) +
  scale_fill_manual("", values = c("Fitted Model" = "black")) +
  xlab("Date") +
  ylab("ILI+ Growth Rate") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
