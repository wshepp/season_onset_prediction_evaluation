# function that calculates ili+
ilip_calc <- function(ili, lci, num_tests){
  return(ili*lci/num_tests)
}

# function that preprocesses raw fluID data using both fluID and FluNet data from WHO
preprocessing <- function(country_name, fluid_df, flunet_df, start_date, end_date){
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
  min_date <- as.Date(start_date)
  max_date <- as.Date(end_date)
  
  fluid_df$week <- as.Date(fluid_df$week)
  flunet_df$week <- as.Date(flunet_df$week)
  
  fluid_df <- fluid_df[fluid_df$week < max_date & fluid_df$week >= min_date, ]
  flunet_df <- flunet_df[flunet_df$week < max_date & flunet_df$week >= min_date, ]
  
  # restricting to country of choice
  fluid_df <- fluid_df[fluid_df$country == country_name, ]
  flunet_df <- flunet_df[flunet_df$country == country_name, ]
  fluid_df$ili[is.na(fluid_df$ili)] <- 0

  # create new dfs removing unused data
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
  
  return(fluid_df)
}


