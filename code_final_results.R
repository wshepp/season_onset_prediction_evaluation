# setup
library(ggplot2)
library(mem)
library(dplyr)
library(mgcv)
library(tidyr)
library(MASS)
EPS <- 0.01
max_weeks <- 442

################################################################################
# Functions
################################################################################

# function that calculates ili+
ilip_calc <- function(ili, lci, num_tests){
  return(ili*lci/num_tests)
}

mem_metrics <- function(data, sth_hem = TRUE){
  metrics_df <- data.frame(
    season = numeric(0),
    epi_season_start = numeric(0),
    peak_week = numeric(0),
    epi_season_end = numeric(0),
    pre_thresh = numeric(0),
    post_thresh = numeric(0),
    false_alerts = numeric(0),
    detection_lag = numeric(0),
    timeliness = numeric(0),
    sens = numeric(0),
    spec = numeric(0),
    ppv = numeric(0),
    npv = numeric(0)
  )
  year_names <- colnames(data)
  for (i in 1:ncol(data)){
    target_season <- data[,year_names[i]]
    hist_seasons <- data[,-which(names(data) == year_names[i])]
    model_curr <- memmodel(i.data = hist_seasons) 
    #i.type.threshold = 6,
    #i.use.t = TRUE)
    pre_thresh_curr <- model_curr$epidemic.thresholds[1]
    season_metrics <- list(season = i,
                           pre_thresh = pre_thresh_curr)
    
    season_df <- data.frame(season = season_metrics$season,
                            pre_thresh = season_metrics$pre_thresh
    )
    metrics_df <- rbind(metrics_df, season_df)
    
  }
  return(metrics_df)
}

# function that finds start of season based on growth rate
find_season_start <- function(season_data, gr_thresh){
  # check for initial negative growth rate
  neg_check <- min(which(season_data < gr_thresh))
  # find positive growth rate after initial wave of negative growth rate
  for (i in neg_check:length(season_data)){
    if (season_data[i] > gr_thresh){
      return(i)
    }
  }
  # no negative growth rate at start? return week 1
  return(1)
}

# function that finds the sason start (kinda useless now and should probs 
# condense into above function)
gr_season_finder <- function(gr_data, start, end, gr_thresh){
  lb_data <- gr_data$lb_90[which(start <= gr_data$time & gr_data$time < end)]
  epi_start <- find_season_start(lb_data, gr_thresh)
  epi_end <- epi_start
  return(list(epi_start = epi_start, epi_end = epi_end))
}

################################################################################
# Start
################################################################################
country_name_vec <- c("Australia", "United States of America")
nth_hem_vec <- c(FALSE, TRUE)
gr_thresh_vec <- c(0, 0.1)
buf_vec <- c(2, 0)
for (country in 1:2){
  country_name <- country_name_vec[country]
  nth_hem <- nth_hem_vec[country]
  gr_thresh <- gr_thresh_vec[country]
  buf <- buf_vec[country]
  
  # read in data
  flunet_df <- read.csv('Data/VIW_FNT.csv', header = T)
  fluid_df <- read.csv('Data/VIW_FID_EPI.csv', header = T)
  
  # renaming columns for ease
  names(fluid_df)[names(fluid_df) == 'ISO_WEEKSTARTDATE'] <- 'week'
  names(flunet_df)[names(flunet_df) == 'ISO_WEEKSTARTDATE'] <- 'week'
  names(fluid_df)[names(fluid_df) == 'ISO_WEEK'] <- 'w_num'
  names(flunet_df)[names(flunet_df) == 'ISO_WEEK'] <- 'w_num'
  names(fluid_df)[names(fluid_df) == 'COUNTRY_AREA_TERRITORY'] <- 'country'
  names(flunet_df)[names(flunet_df) == 'COUNTRY_AREA_TERRITORY'] <- 'country'
  names(fluid_df)[names(fluid_df) == 'ILI_CASE'] <- 'ili'
  names(flunet_df)[names(flunet_df) == 'INF_ALL'] <- 'lci'
  names(flunet_df)[names(flunet_df) == 'INF_A'] <- 'lci_a'
  names(flunet_df)[names(flunet_df) == 'INF_B'] <- 'lci_b'
  names(flunet_df)[names(flunet_df) == 'SPEC_PROCESSED_NB'] <- 'num_tests'
  names(flunet_df)[names(flunet_df) == 'INF_NEGATIVE'] <- 'lci_neg'
  
  # setting time frame
  min_date <- as.Date('2012-01-01')
  max_date <- as.Date('2024-06-30')
  
  fluid_df$week <- as.Date(fluid_df$week)
  flunet_df$week <- as.Date(flunet_df$week)
  
  fluid_df <- fluid_df[fluid_df$week < max_date & fluid_df$week >= min_date, ]
  flunet_df <- flunet_df[flunet_df$week < max_date & flunet_df$week >= min_date, ]
  
  # restricting to country of choice
  fluid_df <- fluid_df[fluid_df$country == country_name, ]
  flunet_df <- flunet_df[flunet_df$country == country_name, ]
  fluid_df$ili[is.na(fluid_df$ili)] <- 0
  
  fluid_df <- summarise(
    group_by(fluid_df, week),
    ili = sum(ili)
  )
  flunet_df <- summarise(
    group_by(flunet_df, week),
    num_tests = sum(num_tests),
    lci = sum(lci),
    lci_a = sum(lci_a),
    lci_b = sum(lci_b),
    lci_neg = sum(lci_neg)
  )
  
  fluid_df <- fluid_df[which(fluid_df$week %in% flunet_df$week),]
  flunet_df <- flunet_df[which(flunet_df$week %in% fluid_df$week),]
  
  # order by date
  fluid_df <- fluid_df[order(fluid_df$week), ]
  flunet_df <- flunet_df[order(flunet_df$week), ]
  
  fluid_df$w_num <- as.integer(strftime(fluid_df$week, format = "%V"))
  flunet_df$w_num <- as.integer(strftime(flunet_df$week, format = "%V"))
  
  # calculating indicators
  
  # lci
  flunet_df$lci_a[is.na(flunet_df$lci_a)] <- 0
  flunet_df$lci_b[is.na(flunet_df$lci_b)] <- 0
  flunet_df$lci[is.na(flunet_df$lci)] <- flunet_df$lci_a[is.na(flunet_df$lci)] + flunet_df$lci_b[is.na(flunet_df$lci)]
  
  # lci negative
  flunet_df$lci_neg[is.na(flunet_df$lci_neg)] <- 0
  
  # tpp
  for (i in which(is.na(flunet_df$num_tests))){
    flunet_df$num_tests[i] <- flunet_df$lci[i] + flunet_df$lci_neg[i]
  }
  flunet_df$tpp <- flunet_df$lci/(flunet_df$num_tests)
  flunet_df$tpp[is.na(flunet_df$tpp)] <- 0
  
  # ili+
  fluid_df$ilip <- fluid_df$ili*flunet_df$tpp
  
  # ili+ for influenza a
  fluid_df$ilip_a <- ilip_calc(fluid_df$ili, flunet_df$lci_a, flunet_df$num_tests)
  # ili+ for influenza b
  fluid_df$ilip_b <- ilip_calc(fluid_df$ili, flunet_df$lci_b, flunet_df$num_tests)
  
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/fluid_df.csv',sep = '')
  write.csv(fluid_df,file_path,row.names = FALSE)
  
  # set starting dates for analysis
  if (nth_hem == TRUE){
    base_date_start <- as.Date("2012-10-01")
    base_date_end <- as.Date("2013-05-27")
  }else{
    base_date_start <- as.Date("2012-03-04")
    base_date_end <- as.Date("2012-11-25")
  }
  
  ################################################################################
  # ILI and ILI+ MEM stuff
  ################################################################################
  
  
  
  # break up ili+ into seasons
  ilip_s1 <- fluid_df$ilip[which(base_date_start <= fluid_df$week & fluid_df$week < base_date_end)]
  ilip_s2 <- fluid_df$ilip[which(base_date_start + 1*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 1*365.25)]
  ilip_s3 <- fluid_df$ilip[which(base_date_start + 2*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 2*365.25)]
  ilip_s4 <- fluid_df$ilip[which(base_date_start + 3*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 3*365.25)]
  ilip_s5 <- fluid_df$ilip[which(base_date_start + 4*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 4*365.25)]
  ilip_s6 <- fluid_df$ilip[which(base_date_start + 5*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 5*365.25)]
  ilip_s7 <- fluid_df$ilip[which(base_date_start + 6*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 6*365.25)]
  ilip_s8 <- fluid_df$ilip[which(base_date_start + 7*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 7*365.25)]
  
  # break up ili into seasons
  ili_s1 <- fluid_df$ili[which(base_date_start <= fluid_df$week & fluid_df$week < base_date_end)]
  ili_s2 <- fluid_df$ili[which(base_date_start + 1*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 1*365.25)]
  ili_s3 <- fluid_df$ili[which(base_date_start + 2*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 2*365.25)]
  ili_s4 <- fluid_df$ili[which(base_date_start + 3*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 3*365.25)]
  ili_s5 <- fluid_df$ili[which(base_date_start + 4*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 4*365.25)]
  ili_s6 <- fluid_df$ili[which(base_date_start + 5*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 5*365.25)]
  ili_s7 <- fluid_df$ili[which(base_date_start + 6*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 6*365.25)]
  ili_s8 <- fluid_df$ili[which(base_date_start + 7*365.25 <= fluid_df$week & fluid_df$week <= base_date_end + 7*365.25)]
  
  # create dataframes
  max_length <- min(length(ilip_s1),
                    length(ilip_s2),
                    length(ilip_s3),
                    length(ilip_s4),
                    length(ilip_s5),
                    length(ilip_s6),
                    length(ilip_s7),
                    length(ilip_s8))
  
  if (nth_hem == FALSE){
    ilip_df <- data.frame('2012' = ilip_s1[1:max_length],
                          '2013' = ilip_s2[1:max_length],
                          '2014' = ilip_s3[1:max_length],
                          '2015' = ilip_s4[1:max_length],
                          '2016' = ilip_s5[1:max_length],
                          '2017' = ilip_s6[1:max_length],
                          '2018' = ilip_s7[1:max_length],
                          '2019' = ilip_s8[1:max_length])
    ili_df <- data.frame('2012' = ili_s1[1:max_length],
                         '2013' = ili_s2[1:max_length],
                         '2014' = ili_s3[1:max_length],
                         '2015' = ili_s4[1:max_length],
                         '2016' = ili_s5[1:max_length],
                         '2017' = ili_s6[1:max_length],
                         '2018' = ili_s7[1:max_length],
                         '2019' = ili_s8[1:max_length])
  }else{
    ilip_df <- data.frame('2012/2013' = ilip_s1[1:max_length],
                          '2013/2014' = ilip_s2[1:max_length],
                          '2014/2015' = ilip_s3[1:max_length],
                          '2015/2016' = ilip_s4[1:max_length],
                          '2016/2017' = ilip_s5[1:max_length],
                          '2017/2018' = ilip_s6[1:max_length],
                          '2018/2019' = ilip_s7[1:max_length],
                          '2019/2020' = ilip_s8[1:max_length])
    ili_df <- data.frame('2012/2013' = ili_s1[1:max_length],
                         '2013/2014' = ili_s2[1:max_length],
                         '2014/2015' = ili_s3[1:max_length],
                         '2015/2016' = ili_s4[1:max_length],
                         '2016/2017' = ili_s5[1:max_length],
                         '2017/2018' = ili_s6[1:max_length],
                         '2018/2019' = ili_s7[1:max_length],
                         '2019/2020' = ili_s8[1:max_length])
  }
  
  # save data
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/ili_df.csv',sep = '')
  write.csv(ili_df,file_path,row.names = FALSE)
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/ilip_df.csv',sep = '')
  write.csv(ilip_df,file_path,row.names = FALSE)
  
  # calculate and save metrics
  ili_metrics <- mem_metrics(ili_df, sth_hem = TRUE)
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/ili_thresholds.csv',sep = '')
  write.csv(ili_metrics,file_path,row.names = FALSE)
  ilip_metrics <- mem_metrics(ilip_df, sth_hem = TRUE)
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/ili_plus_thresholds.csv',sep = '')
  write.csv(ilip_metrics,file_path,row.names = FALSE)
  
  ################################################################################
  # Data and Graph Generating - Run time ~10 min
  ################################################################################
  
  set.seed(123)
  for (i in 5:max_weeks){
    # read in data for current plot
    data_df <- data.frame(
      week = fluid_df$week[1:i],
      ilip = fluid_df$ilip[1:i]
    )
    # generate model for each indicator
    model <- gam(log(ilip+EPS)~s(as.numeric(week), k = as.integer(floor(i/5))), 
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
    for (j in 1:i){
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
    # save current data
    file_path <- file_path <- paste('Final Results/GAM/',country_name,'/data/cases/week_',i,'.csv',sep = '')
    write.csv(plotting_df,file_path,row.names = FALSE)
    
    # plot current graph
    ggplot(data = plotting_df) +
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
      geom_point(data = data_df, aes(x = week, y = ilip)) +
      coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                               as.Date("2020-07-01") - 180)) +
      scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
      scale_color_manual("", values = c("Fitted Model" = "black")) +
      scale_fill_manual("", values = c("Fitted Model" = "black")) +
      xlab("Date") +
      ylab("ILI+ Values") +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")
    file_path <- paste('Final Results/GAM/',country_name,'/graphs/cases/week_',i,'.png',sep = '')
    ggsave(file_path, width = 6.67, height = 5, units = "in")
    
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
    for (j in 1:i){
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
    
    # save current data
    file_path <- file_path <- paste('Final Results/GAM/',country_name,'/data/growth_rate/week_',i,'.csv',sep = '')
    write.csv(plotting_gr_df,file_path,row.names = FALSE)
    
    # plot current graph
    ggplot(data = plotting_gr_df) +
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
      geom_hline(yintercept = gr_thresh,
                 color = 'red4',
                 linetype = 2) +
      coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                               as.Date("2020-07-01") - 180)) +
      scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
      scale_color_manual("", values = c("Fitted Model" = "black")) +
      scale_fill_manual("", values = c("Fitted Model" = "black")) +
      xlab("Date") +
      ylab("ILI+ Growth Rate") +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")
    file_path <- paste('Final Results/GAM/',country_name,'/graphs/growth_rate/week_',i,'.png',sep = '')
    ggsave(file_path, width = 6.67, height = 5, units = "in")
  }
  
  ################################################################################
  
  # generate real time growth rate data
  rt_gr <- data.frame()
  for (i in 5:max_weeks){
    file_path <- paste('Final Results/GAM/',country_name,'/data/growth_rate/week_',i,'.csv',sep = '')
    curr_model <- read.csv(file_path, header = TRUE)
    curr_model$time <- as.Date(curr_model$time)
    data_length <- length(curr_model$time)
    row <- data.frame(time = curr_model$time[data_length],
                      y = curr_model$y[data_length],
                      lb_50 = curr_model$lb_50[data_length],
                      ub_50 = curr_model$ub_50[data_length],
                      lb_80 = curr_model$lb_80[data_length],
                      ub_80 = curr_model$ub_80[data_length],
                      lb_90 = curr_model$lb_80[data_length],
                      ub_90 = curr_model$ub_80[data_length],
                      lb_95 = curr_model$lb_95[data_length],
                      ub_95 = curr_model$ub_95[data_length])
    rt_gr <- rbind(rt_gr, row)
  }
  file_path <- paste('Final Results/Cleaned Data/',country_name,'/rt_gr_model.csv',sep = '')
  write.csv(rt_gr,file_path,row.names = FALSE)
  
  # set year labels based on hemisphere
  if (nth_hem == TRUE){
    year_labels <- c(rep("2012/2013", 2),
                     rep("2013/2014", 2),
                     rep("2014/2015", 2),
                     rep("2015/2016", 2),
                     rep("2016/2017", 2),
                     rep("2017/2018", 2),
                     rep("2018/2019", 2),
                     rep("2019/2020", 2))
  }else{
    year_labels <- c(rep("2012", 2),
                     rep("2013", 2),
                     rep("2014", 2),
                     rep("2015", 2),
                     rep("2016", 2),
                     rep("2017", 2),
                     rep("2018", 2),
                     rep("2019", 2))
  }
  
  # create base weeks df to populate
  weeks_data <- data.frame(
    indicator = rep(c("ILI+ Quans",
                      "ILI+ GAM GR - RT"),8),
    year = year_labels,
    start = rep(as.Date("2012-01-01"),16),
    end = rep(as.Date("2012-01-01"),16)
  )
  
  ################################################################################
  # find starts to seasons as well as coverage area
  for (i in 1:8){
    s_curr_start <- base_date_start + (i - 1)*365.25
    s_curr_end <- base_date_end + (i - 1)*365.25
    start <- fluid_df$week[which(s_curr_start < fluid_df$week)][1]
    end <- tail(fluid_df$week[which(s_curr_end >= fluid_df$week)], n=1)
    ilip_curr <- fluid_df$ilip[which(start <= fluid_df$week & fluid_df$week < end)]
    csum_curr <- cumsum(ilip_curr[5:length(ilip_curr)]/sum(ilip_curr[5:length(ilip_curr)]))
    csum_start_curr <-min(which(csum_curr > 0.025))+4
    csum_end_curr <- min(which(csum_curr > 0.975))+4
    if (nth_hem == TRUE){
      csum_curr <- cumsum(ilip_curr/sum(ilip_curr))
      csum_start_curr <-min(which(csum_curr > 0.025))
      csum_end_curr <- min(which(csum_curr > 0.975))
    }
    gr_season_info_rt_curr <- gr_season_finder(rt_gr, s_curr_start+buf*7, s_curr_end, gr_thresh)
    s_start_curr <- c(csum_start_curr,
                      gr_season_info_rt_curr$epi_start+buf)
    s_end_curr <- c(csum_end_curr,
                    gr_season_info_rt_curr$epi_start+buf)
    for (j in 1:2){
      weeks_data$start[(i-1)*2+j] <- as.Date(start+7*(s_start_curr[j]-1)-(i-1)*365.25)
      weeks_data$end[(i-1)*2+j] <- as.Date(start+7*(s_end_curr[j]-1)-(i-1)*365.25)
    }
  }
  
  # plot season graphs
  pred_eval_df <- data.frame()
  pred_data <- data.frame()
  if (country_name == "Australia") {
    for (sn in 1:8){
      s_curr_start <- base_date_start + (sn - 1)*365.25
      s_curr_end <- base_date_end + (sn - 1)*365.25
      s_week_start <- fluid_df$week[which(s_curr_start < fluid_df$week)][1]
      s_week_end <- s_week_start + 7*(max_length-1)
      s_curr <- data.frame(
        week = seq(s_week_start,s_week_end,7),
        ili = ili_df[,sn],
        ili_plus = ilip_df[,sn]
      )
      quan_dates <- data.frame(start = (weeks_data$start[(sn-1)*2+1]+(sn-1)*365.25),
                               end = (weeks_data$end[(sn-1)*2+1]+(sn-1)*365.25))
      rt_gr_start <- data.frame(start = (weeks_data$start[(sn-1)*2+2]+(sn-1)*365.25))
      mem_pred_start_ili <- s_curr$week[min(which(s_curr$ili > ili_metrics$pre_thresh[s_num]))]
      mem_pred_start_ilip <- s_curr$week[min(which(s_curr$ili_plus > ilip_metrics$pre_thresh[s_num]))]
      ili_start <- data.frame(start = mem_pred_start_ili)
      ilip_start <- data.frame(start = mem_pred_start_ilip)
      
      # season start prediction comparion plot
      ggplot() +
        geom_vline(aes(xintercept = start, color = "MEM - ILI"),
                   linetype = 2,
                   linewidth = 2,
                   data = ili_start) +
        geom_vline(aes(xintercept = start, color = "MEM - ILI+"),
                   linetype = 3,
                   linewidth = 4,
                   data = ilip_start) +
        geom_rect(aes(xmin = start,xmax = end,
                      ymin = -Inf,ymax = Inf,
                      fill = "Central 95% Coverage"),
                  alpha = 0.2,
                  data = quan_dates) +
        geom_vline(aes(xintercept = start, color = "GAM Growth Rate"),
                   linetype = 6,
                   linewidth = 2,
                   data = rt_gr_start) +
        geom_line(aes(x = week, y = ili_plus), 
                  color = "navy", 
                  linewidth = 2,
                  data = s_curr) +
        geom_point(aes(x = week, y = ili_plus),
                   fill = NA,
                   color = "black",
                   size = 5,
                   data = s_curr) +
        
        labs(x = "Date",
             y = "ILI+ Value") +
        theme_bw(base_size = 39) +
        theme(legend.position = "bottom") + 
        scale_fill_manual(
          name = "",
          values = c("Central 95% Coverage" = "grey40")) + 
        scale_color_manual(
          name = "",
          values = c("MEM - ILI" = "navy",
                     "MEM - ILI+" = "green4",
                     "GAM Growth Rate" = "purple3")) +
        scale_x_date(
          limits = c(s_week_start,s_week_end),
          date_breaks = "1 month",
          date_labels = "%b"
        )
      file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/season_',sn,'.png',sep = '')
      ggsave(file_path, width = 18, height = 12, units = "in")
      
      # calculate eval metrics for each season
      ilip_curr <- fluid_df$ilip[which(s_week_start <= fluid_df$week & fluid_df$week < s_week_end)]
      csum_curr <- cumsum(ilip_curr[5:length(ilip_curr)]/sum(ilip_curr[5:length(ilip_curr)]))
      start_week_n <- (min(which(csum_curr > 0.025)) + 4)
      if(nth_hem == TRUE){
        csum_curr <- cumsum(ilip_curr/sum(ilip_curr))
        start_week_n <- (min(which(csum_curr > 0.025)))
      }
      ili_weeks_late <- as.numeric(ili_start$start - quan_dates$start)/7
      ilip_weeks_late <- as.numeric(ilip_start$start - quan_dates$start)/7
      gr_weeks_late <- as.numeric(rt_gr_start$start - quan_dates$start)/7
      if (ili_weeks_late < 0 | is.na(ili_weeks_late)){
        ili_coverage_missed <- NA
      }else{
        ili_coverage_missed <- (csum_curr[start_week_n+ili_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      if (ilip_weeks_late < 0 | is.na(ilip_weeks_late)){
        ilip_coverage_missed <- NA
      }else{
        ilip_coverage_missed <- (csum_curr[start_week_n+ilip_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      if (gr_weeks_late < 0 | is.na(gr_weeks_late)){
        gr_coverage_missed <- NA
      }else{
        gr_coverage_missed <- (csum_curr[start_week_n+gr_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      
      row <- data.frame(
        season = sn,
        start_week = quan_dates$start,
        end_week = quan_dates$end,
        start_week_n = start_week_n,
        ili_weeks_late = as.numeric(mem_pred_start_ili - quan_dates$start)/7,
        ili_coverage_missed = round(ili_coverage_missed,3),
        ilip_weeks_late = as.numeric(mem_pred_start_ilip - quan_dates$start)/7,
        ilip_coverage_missed = round(ilip_coverage_missed,3),
        gr_weeks_late = as.numeric(rt_gr_start$start - quan_dates$start)/7,
        gr_coverage_missed = round(gr_coverage_missed,3)
      )
      pred_eval_df <- rbind(pred_eval_df, row)
      
      row <- data.frame(
        year = weeks_data$year[(sn-1)*2+1],
        coverage_start = weeks_data$start[(sn-1)*2+1],
        coverage_end = weeks_data$end[(sn-1)*2+1],
        ili_start = mem_pred_start_ili-365.25*(sn-1),
        ilip_start = mem_pred_start_ilip-365.25*(sn-1),
        gr_rt_start = weeks_data$start[(sn-1)*2+2]
      )
      pred_data <- rbind(pred_data, row)
    }
    
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/prediction_eval_metrics.csv',sep = '')
    write.csv(pred_eval_df,file_path,row.names = FALSE)
    
    # summarise eval metrics
    eval_tot_df <- data.frame(
      predictor = c("MEM - ILI", 
                    "MEM - ILI+", 
                    "GAM Growth Rate - ILI+"),
      avg_coverage_missed = c(round(mean(na.omit(pred_eval_df$ili_coverage_missed)),3),
                              round(mean(na.omit(pred_eval_df$ilip_coverage_missed)),3),
                              round(mean(na.omit(pred_eval_df$gr_coverage_missed)),3)),
      median_late_weeks = c(median(pred_eval_df$ili_weeks_late[which(pred_eval_df$ili_weeks_late >= 0)]),
                            median(pred_eval_df$ilip_weeks_late[which(pred_eval_df$ilip_weeks_late >= 0)]),
                            median(pred_eval_df$gr_weeks_late[which(pred_eval_df$gr_weeks_late >= 0)])),
      median_early_weeks = c(-median(pred_eval_df$ili_weeks_late[which(pred_eval_df$ili_weeks_late <= 0)]),
                             -median(pred_eval_df$ilip_weeks_late[which(pred_eval_df$ilip_weeks_late <= 0)]),
                             -median(pred_eval_df$gr_weeks_late[which(pred_eval_df$gr_weeks_late <= 0)])),
      seasons_missed = c(sum(is.na(pred_eval_df$ili_weeks_late)),
                         sum(is.na(pred_eval_df$ilip_weeks_late)),
                         sum(is.na(pred_eval_df$gr_weeks_late)))
    )
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/prediction_eval_summary.csv',sep = '')
    write.csv(eval_tot_df,file_path,row.names = FALSE)
    
    ggplot(pred_data, aes(x = coverage_start, xend = coverage_end, y = year, yend = year)) +
      geom_segment(aes(color = "Central 95% Coverage"), size = 1, alpha = 0.5) + 
      geom_point(aes(x = ili_start, y = year, color = "ILI MEM Predicted Start"),
                 size = 3,
                 shape = 8,
                 stroke = 1) +
      geom_point(aes(x = ilip_start, y = year, color = "ILI+ MEM Predicted Start"),
                 size = 3,
                 shape = 4,
                 stroke = 1) +
      geom_point(aes(x = gr_rt_start, y = year, color = "GAM Real Time Predicted Start"),
                 size = 3,
                 shape = 1,
                 stroke = 1) +
      labs(title = "Onset Predictions By Year", x = "Date", y = "Year") + 
      scale_x_date(
        limits = c(base_date_start,base_date_end),
        date_breaks = "1 month",
        date_labels = "%b"
      ) + 
      scale_color_manual(
        name = "",
        values = c("ILI MEM Predicted Start" = "navy",
                   "ILI+ MEM Predicted Start" = "green4",
                   "GAM Real Time Predicted Start" = "purple3",
                   "Central 95% Coverage" = "red4") 
      ) +
      theme_minimal(base_size = 22) +
      theme(legend.position = "bottom")
    
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/predictions_summary.png',sep = '')
    ggsave(file_path, width = 10, height = 10, units = "in")
  }
  
   ################################################################################
  
  if (country_name == "United States of America") {
    usa_dates <- c(as.Date("2012-12-17"),
                   as.Date("2013-12-09"),
                   as.Date("2014-11-24"),
                   as.Date("2015-12-21"),
                   as.Date("2016-12-19"),
                   as.Date("2017-11-27"),
                   as.Date("2018-12-17"),
                   as.Date("2019-11-25"))
    
    for (sn in 1:8){
      s_curr_start <- base_date_start + (sn - 1)*365.25
      s_curr_end <- base_date_end + (sn - 1)*365.25
      s_week_start <- fluid_df$week[which(s_curr_start < fluid_df$week)][1]
      s_week_end <- s_week_start + 7*(max_length-1)
      s_curr <- data.frame(
        week = seq(s_week_start,s_week_end,7),
        ili = ili_df[,sn],
        ili_plus = ilip_df[,sn]
      )
      quan_dates <- data.frame(start = (weeks_data$start[(sn-1)*2+1]+(sn-1)*365.25),
                               end = (weeks_data$end[(sn-1)*2+1]+(sn-1)*365.25))
      rt_gr_start <- data.frame(start = (weeks_data$start[(sn-1)*2+2]+(sn-1)*365.25))
      mem_pred_start_ili <- s_curr$week[min(which(s_curr$ili > ili_metrics$pre_thresh[s_num]))]
      mem_pred_start_ilip <- s_curr$week[min(which(s_curr$ili_plus > ilip_metrics$pre_thresh[s_num]))]
      ili_start <- data.frame(start = mem_pred_start_ili)
      ilip_start <- data.frame(start = mem_pred_start_ilip)
      usa_start <- data.frame(start = usa_dates[sn])
      
      # season start prediction comparion plot
      ggplot() +
        geom_vline(aes(xintercept = start, color = "MEM - ILI"),
                   linetype = 2,
                   linewidth = 0.75,
                   data = ili_start) +
        geom_vline(aes(xintercept = start, color = "MEM - ILI+"),
                   linetype = 3,
                   linewidth = 0.75,
                   data = ilip_start) +
        geom_vline(aes(xintercept = start, color = "USA ILI%"),
                   linetype = 4,
                   linewidth = 0.75,
                   data = usa_start) +
        geom_rect(aes(xmin = start,xmax = end,
                      ymin = -Inf,ymax = Inf,
                      fill = "Central 95% Coverage"),
                  alpha = 0.2,
                  data = quan_dates) +
        geom_vline(aes(xintercept = start, color = "GAM Growth Rate"),
                   linetype = 6,
                   linewidth = 0.75,
                   data = rt_gr_start) +
        geom_line(aes(x = week, y = ili_plus), 
                  color = "navy", 
                  data = s_curr) +
        geom_point(aes(x = week, y = ili_plus),
                   fill = NA,
                   color = "black",
                   size = 2,
                   data = s_curr) +
        
        labs(x = "Date",
             y = "ILI+ Value") +
        theme_bw(base_size = 22) +
        theme(legend.position = "bottom") + 
        scale_fill_manual(
          name = "",
          values = c("Central 95% Coverage" = "grey40")) + 
        scale_color_manual(
          name = "",
          values = c("MEM - ILI" = "navy",
                     "MEM - ILI+" = "green4",
                     "USA ILI%" = "red4",
                     "GAM Growth Rate" = "purple3")) +
        scale_x_date(
          limits = c(s_week_start,s_week_end),
          date_breaks = "1 month",
          date_labels = "%b"
        )
      file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/season_',sn,'.png',sep = '')
      ggsave(file_path, width = 13.33, height = 10, units = "in")
      
      ilip_curr <- fluid_df$ilip[which(s_week_start <= fluid_df$week & fluid_df$week < s_week_end)]
      csum_curr <- cumsum(ilip_curr[5:length(ilip_curr)]/sum(ilip_curr[5:length(ilip_curr)]))
      start_week_n <- (min(which(csum_curr > 0.025)) + 4)
      if(nth_hem == TRUE){
        csum_curr <- cumsum(ilip_curr/sum(ilip_curr))
        start_week_n <- (min(which(csum_curr > 0.025)))
      }
      ili_weeks_late <- as.numeric(ili_start$start - quan_dates$start)/7
      ilip_weeks_late <- as.numeric(ilip_start$start - quan_dates$start)/7
      usa_weeks_late <- as.numeric(usa_dates[sn] - quan_dates$start)/7
      gr_weeks_late <- as.numeric(rt_gr_start$start - quan_dates$start)/7
      if (ili_weeks_late < 0 | is.na(ili_weeks_late)){
        ili_coverage_missed <- NA
      }else{
        ili_coverage_missed <- (csum_curr[start_week_n+ili_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      if (ilip_weeks_late < 0 | is.na(ilip_weeks_late)){
        ilip_coverage_missed <- NA
      }else{
        ilip_coverage_missed <- (csum_curr[start_week_n+ilip_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      if (usa_weeks_late < 0 | is.na(usa_weeks_late)){
        usa_coverage_missed <- NA
      }else{
        usa_coverage_missed <- (csum_curr[start_week_n+usa_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      if (gr_weeks_late < 0 | is.na(gr_weeks_late)){
        gr_coverage_missed <- NA
      }else{
        gr_coverage_missed <- (csum_curr[start_week_n+gr_weeks_late] - csum_curr[start_week_n])/0.95*100
      }
      
      row <- data.frame(
        season = sn,
        start_week = quan_dates$start,
        end_week = quan_dates$end,
        start_week_n = start_week_n,
        ili_weeks_late = as.numeric(mem_pred_start_ili - quan_dates$start)/7,
        ili_coverage_missed = round(ili_coverage_missed,3),
        ilip_weeks_late = as.numeric(mem_pred_start_ilip - quan_dates$start)/7,
        ilip_coverage_missed = round(ilip_coverage_missed,3),
        usa_weeks_late = usa_weeks_late,
        usa_coverage_missed = round(usa_coverage_missed,3),
        gr_weeks_late = as.numeric(rt_gr_start$start - quan_dates$start)/7,
        gr_coverage_missed = round(gr_coverage_missed,3)
      )
      pred_eval_df <- rbind(pred_eval_df, row)
      
      row <- data.frame(
        year = weeks_data$year[(sn-1)*2+1],
        coverage_start = weeks_data$start[(sn-1)*2+1],
        coverage_end = weeks_data$end[(sn-1)*2+1],
        ili_start = mem_pred_start_ili-365.25*(sn-1),
        ilip_start = mem_pred_start_ilip-365.25*(sn-1),
        usa_start = usa_dates[sn]-365.25*(sn-1),
        gr_rt_start = weeks_data$start[(sn-1)*2+2]
      )
      pred_data <- rbind(pred_data, row)
    }
    
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/prediction_eval_metrics.csv',sep = '')
    write.csv(pred_eval_df,file_path,row.names = FALSE)
    eval_tot_df <- data.frame(
      predictor = c("MEM - ILI", 
                    "MEM - ILI+",
                    "USA ILI%",
                    "GAM Growth Rate - ILI+"),
      avg_coverage_missed = c(round(mean(na.omit(pred_eval_df$ili_coverage_missed)),3),
                              round(mean(na.omit(pred_eval_df$ilip_coverage_missed)),3),
                              round(mean(na.omit(pred_eval_df$usa_coverage_missed)),3),
                              round(mean(na.omit(pred_eval_df$gr_coverage_missed)),3)),
      median_late_weeks = c(median(pred_eval_df$ili_weeks_late[which(pred_eval_df$ili_weeks_late >= 0)]),
                            median(pred_eval_df$ilip_weeks_late[which(pred_eval_df$ilip_weeks_late >= 0)]),
                            median(pred_eval_df$usa_weeks_late[which(pred_eval_df$usa_weeks_late >= 0)]),
                            median(pred_eval_df$gr_weeks_late[which(pred_eval_df$gr_weeks_late >= 0)])),
      median_early_weeks = c(-median(pred_eval_df$ili_weeks_late[which(pred_eval_df$ili_weeks_late <= 0)]),
                             -median(pred_eval_df$ilip_weeks_late[which(pred_eval_df$ilip_weeks_late <= 0)]),
                             -median(pred_eval_df$usa_weeks_late[which(pred_eval_df$usa_weeks_late <= 0)]),
                             -median(pred_eval_df$gr_weeks_late[which(pred_eval_df$gr_weeks_late <= 0)])),
      seasons_missed = c(sum(is.na(pred_eval_df$ili_weeks_late)),
                         sum(is.na(pred_eval_df$ilip_weeks_late)),
                         sum(is.na(pred_eval_df$usa_weeks_late)),
                         sum(is.na(pred_eval_df$gr_weeks_late)))
    )
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/prediction_eval_summary.csv',sep = '')
    write.csv(eval_tot_df,file_path,row.names = FALSE)
    
    ggplot(pred_data, aes(x = coverage_start, xend = coverage_end, y = year, yend = year)) +
      geom_segment(aes(color = "Central 95% Coverage"), size = 1, alpha = 0.5) + 
      geom_point(aes(x = ili_start, y = year, color = "ILI MEM Predicted Start"),
                 size = 3,
                 shape = 8,
                 stroke = 1) +
      geom_point(aes(x = ilip_start, y = year, color = "ILI+ MEM Predicted Start"),
                 size = 3,
                 shape = 4,
                 stroke = 1) +
      geom_point(aes(x = usa_start, y = year, color = "USA ILI% Predicted Start"),
                 size = 3,
                 shape = 5,
                 stroke = 1) +
      geom_point(aes(x = gr_rt_start, y = year, color = "GAM Real Time Predicted Start"),
                 size = 3,
                 shape = 1,
                 stroke = 1) +
      labs(title = "Onset Predictions By Year", x = "Date", y = "Year") + 
      scale_x_date(
        limits = c(base_date_start,base_date_end),
        date_breaks = "1 month",
        date_labels = "%b"
      ) + 
      scale_color_manual(
        name = "",
        values = c("ILI MEM Predicted Start" = "navy",
                   "ILI+ MEM Predicted Start" = "green4",
                   "GAM Real Time Predicted Start" = "purple3",
                   "Central 95% Coverage" = "red4") 
      ) +
      theme_minimal(base_size = 22) +
      theme(legend.position = "bottom")
    
    file_path <- paste('Final Results/',country_name,'/Predictive Season Starts Comparisons/predictions_summary.png',sep = '')
    ggsave(file_path, width = 10, height = 10, units = "in")
  } 
}
################################################################################

# combine eval datasets
pred_eval_df_aus <- read.csv('Final Results/Australia/Predictive Season Starts Comparisons/prediction_eval_metrics.csv', header = T)
pred_eval_df_usa <- read.csv('Final Results/United States of America/Predictive Season Starts Comparisons/prediction_eval_metrics.csv', header = T)
weeks_df_aus <- data.frame(
  predictor = c(rep("MEM - ILI",8), 
                rep("MEM - ILI+",8),
                rep("GAM Growth Rate",8)),
  country = c(rep("Australia", 24)),
  weeks_late = c(pred_eval_df_aus$ili_weeks_late,
                 pred_eval_df_aus$ilip_weeks_late,
                 pred_eval_df_aus$gr_weeks_late)
)
weeks_df_usa <- data.frame(
  predictor = c(rep("MEM - ILI",8), 
                rep("MEM - ILI+",8),
                rep("USA ILI%",8),
                rep("GAM Growth Rate",8)),
  country = c(rep("United States of America", 32)),
  weeks_late = c(pred_eval_df_usa$ili_weeks_late,
                 pred_eval_df_usa$ilip_weeks_late,
                 pred_eval_df_usa$usa_weeks_late,
                 pred_eval_df_usa$gr_weeks_late)
)
weeks_df <- rbind(weeks_df_aus, weeks_df_usa)

cols <- c("navy", "red4" )
weeks_df$predictor <- as.factor(weeks_df$predictor)


median_weeks <- data.frame(
  predictor = c("MEM - ILI", 
                "MEM - ILI+",
                "USA ILI%",
                "GAM Growth Rate"),
  median_weeks_late = c(median(weeks_df$weeks_late[which(weeks_df$predictor == "MEM - ILI")]),
                        median(weeks_df$weeks_late[which(weeks_df$predictor == "MEM - ILI+")]),
                        median(weeks_df$weeks_late[which(weeks_df$predictor == "USA ILI%")]),
                        median(weeks_df$weeks_late[which(weeks_df$predictor == "GAM Growth Rate")])),
  average_weeks_late = c(mean(weeks_df$weeks_late[which(weeks_df$predictor == "MEM - ILI")]),
                         mean(weeks_df$weeks_late[which(weeks_df$predictor == "MEM - ILI+")]),
                         mean(weeks_df$weeks_late[which(weeks_df$predictor == "USA ILI%")]),
                         mean(weeks_df$weeks_late[which(weeks_df$predictor == "GAM Growth Rate")]))
)
# plot prediction distribution
set.seed(123)
ggplot(weeks_df, aes(y = predictor)) +
  geom_boxplot(
    aes(x = weeks_late),
    outlier.shape = NA,
    alpha = 0.4,
    size = 1.3,
    na.rm = TRUE
  ) +
  geom_jitter(
    data = subset(weeks_df, !is.na(weeks_late)),
    aes(x = weeks_late, color = country),
    width = 0.6,
    height = 0,
    alpha = 0.4,
    size = 4
  ) +
  geom_jitter(
    data = subset(weeks_df, is.na(weeks_late)),
    aes(x = 16, color = country),
    width = 0.6,
    height = 0,
    alpha = 0.4,
    size = 4
  ) +
  geom_vline(xintercept = 0) +
  labs(x = "Weeks from Season Start", y = "Method") +
  scale_color_manual(name = "", values = cols) +
  scale_x_continuous(
    breaks = c(-6,-4,-2,0,2,4,6,8,10,12,16),
    labels = c(-6,-4,-2,0,2,4,6,8,10,12,"Missed")
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.y = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
file_path <- 'Final Results/Final Figures/weeks_late_summary.png'
ggsave(file_path, width = 10, height = 8, units = "in")

covarge_missed_df <- data.frame(
  predictor = c("MEM - ILI", 
                "MEM - ILI+", 
                "GAM Growth Rate",
                "MEM - ILI", 
                "MEM - ILI+", 
                "USA ILI%",
                "GAM Growth Rate"),
  country = c(rep("Australia",3),
              rep("United States of America", 4)),
  avg_coverage_missed = c(round(mean(na.omit(pred_eval_df_aus$ili_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_aus$ilip_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_aus$gr_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_usa$ili_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_usa$ilip_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_usa$usa_coverage_missed)),3),
                          round(mean(na.omit(pred_eval_df_usa$gr_coverage_missed)),3))
)

cols = c("navy", "red4")

# plot average coverage missed
ggplot(covarge_missed_df, aes(x = predictor, y = avg_coverage_missed, fill = country)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(),
    alpha = 0.8,
    color = "grey30",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = paste0(round(avg_coverage_missed,1),'%')),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 13.4
  ) +
  scale_fill_manual(name = "", values = cols) + 
  theme_minimal(base_size = 56) +
  labs(x = "Method", y = "Average Coverage Missed") +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

file_path <- 'Final Results/Final Figures/coverage_missed_summary.png'
ggsave(file_path, width = 20, height = 16, units = "in")

