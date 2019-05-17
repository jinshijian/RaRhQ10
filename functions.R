#*****************************************************************************************************************
# Baasic functions
#*****************************************************************************************************************
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=F,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     #do.call("rbind", tapply(xx[[col]], measurevar, quantile, c(0.25, 0.5, 0.75)))
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# filtration function
filtration <- function (sdata, Q10_type, RC) {
  n_col <- colnames(sdata) == Q10_type
  n_RC <- colnames(sdata) == RC
  sdata <- sdata[!is.na(sdata[, n_col]), ]
  sdata <- sdata[!is.na(sdata[, n_RC]) & sdata[, n_RC] > 0, ]
  sdata <- subset(sdata, select=c("Record_number", "Study_number", "Species", "Biome", "Ecosystem_type", "Leaf_habit",
                                  "R10", Q10_type, RC  ) )
  sdata$Q10_type <- Q10_type
  colnames(sdata) <- c( "Record_number", "Study_number", "Species", "Biome", "Ecosystem_type", "Leaf_habit",
                        "R10", "Q10", "RC_annual", "Q10_type" )
  
  return(sdata)
}

#*****************************************************************************************************************
# functions for SRDB 
#*****************************************************************************************************************
# test whether Q10 differ from Ra and Rh dominated sites
Q10_test <- function(sdata, Q10_type, var_title) {
  print(SEPARATOR)
  n_col <- colnames(sdata) == Q10_type
  sdata <- sdata[!is.na(sdata[,n_col]), ]
  sdata <- sdata[!is.na(sdata$RC_annual), ]
  sdata$Q10 <- sdata[, colnames(sdata) == Q10_type]
  sdata <- sdata[sdata$Q10 > 0, ]
  sdata$Ra_Rh_dom <- ifelse(sdata$RC_annual>0.5, 'Ra_dom', 'Rh_dom')
  n_Ra <- nrow(subset(sdata, sdata$Ra_Rh_dom=='Ra_dom'))
  n_Rh <- nrow(subset(sdata, sdata$Ra_Rh_dom=='Rh_dom'))
  
  print(paste0('**********normal-test**********', Q10_type))
  print(shapiro.test(sdata[sdata$Ra_Rh_dom == "Ra_dom", ]$Q10)) # test normality of Ra
  print(shapiro.test(sdata[sdata$Ra_Rh_dom == "Rh_dom", ]$Q10)) # test normality of Rh
  # test whether varance differ
  print(var.test(Q10 ~ Ra_Rh_dom, data = sdata))
  
  print(paste0('**********ANOVA**********'))
  anova <- aov(Q10 ~ Ra_Rh_dom, data = sdata)
  print (summary(anova))
  
  print(paste0('**********t-test**********'))
  t_t <- t.test(Q10 ~ Ra_Rh_dom, data = sdata, var.eq=TRUE)
  print( t_t ) 
  
  print(paste0('**********wilcox-test**********'))
  # boxplot(sdata$Q10 ~ sdata$Ra_Rh_dom)
  wilcox <- wilcox.test(Q10 ~ Ra_Rh_dom, data = sdata, exact=F)
  print( wilcox ) 
  
  p <- ggplot(sdata, aes(Q10, color=RC_annual>0.5, fill = RC_annual>0.5) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(var_title) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    scale_fill_discrete(name="Ra dominate sites",labels = c(paste0("FALSE (n=", n_Ra, ")"), paste0("TRUE (n=", n_Rh, ")") ) ) +
    scale_color_discrete(name="Ra dominate sites",labels = c(paste0("FALSE (n=", n_Ra, ")"), paste0("TRUE (n=", n_Rh, ")") ) )
}

R10_test <- function (sdata, R10_type, var_title) {
  print(SEPARATOR)
  sdata <- sdata[!is.na(sdata$RC_annual), ]
  sdata <- sdata[!is.na(sdata$R10), ]
  sdata$R10 <- sdata[, colnames(sdata) == R10_type]
  sdata$Ra_Rh_dom <- ifelse(sdata$RC_annual>0.5, 'Ra_dom', 'Rh_dom')
  n_Ra <- nrow(subset(sdata, sdata$Ra_Rh_dom=='Ra_dom'))
  n_Rh <- nrow(subset(sdata, sdata$Ra_Rh_dom=='Rh_dom'))
  
  p <- ggplot(sdata, aes(R10, color=RC_annual>0.5, fill = RC_annual>0.5) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(var_title) + xlab(expression( R[10] )) +
    scale_fill_discrete(name="Ra dominate sites",labels = c(paste0("FALSE (n=", n_Ra, ")"), paste0("TRUE (n=", n_Rh, ")") ) )+
    scale_color_discrete(name="Ra dominate sites",labels = c(paste0("FALSE (n=", n_Ra, ")"), paste0("TRUE (n=", n_Rh, ")") ) )
  print(p)
}


#*********************************************************************************
# test Q10 change with SPI and pdsi

Q10_RC_drought <- function (sdata, Q10_type, drought, RC) {
  n_col <- colnames(sdata) == Q10_type
  sdata <- sdata[!is.na(sdata[, n_col]), ]
  sdata <- sdata[sdata[, n_col] < 10, ]
  sdata <- sdata[sdata$Leaf_habit == 'Deciduous' | sdata$Leaf_habit == 'Evergreen', ]
  sdata <- sdata[sdata$Ecosystem_type == 'Forest' | sdata$Ecosystem_type == 'Grassland'
                 | sdata$Ecosystem_type == 'Agriculture' | sdata$Ecosystem_type == 'Savanna'
                 | sdata$Ecosystem_type == 'Woodland', ]
  
  n_RC <- colnames(sdata) == RC
  # sdata <- sdata[!is.na(sdata[, n_RC]) & sdata[, n_RC] > 0, ]
  
  sdata <- subset(sdata, select=c("Record_number", "Study_number", "Species", "Biome", "Ecosystem_type", "Leaf_habit",
                                  "R10", Q10_type, drought, RC  ) )
  
  sdata$Drought_label <- drought
  sdata$Ra_Rh_dom <- ifelse(sdata$RC_annual>0.5, 'Ra_dom', 'Rh_dom')
  colnames(sdata) <- c( "Record_number", "Study_number", "Species", "Biome", "Ecosystem_type", "Leaf_habit",
                        "R10", "Q10",'Drought', "RC_annual", "Drought_label", "Ra_Rh_dom" )
  
  return(sdata)
}



#*****************************************************************************************************************
# functions for MGRsD
#*****************************************************************************************************************
# modeling and calculate Q10
# colnames(MGRsD)
RhRa_M <- function (sdata) {
  sub_Rh <- subset(sdata, !is.na(Rh), select = c(StudyNumber, Rh, Tm ))
  sub_Rh$Source <- "MGRsD_TA"
  sub_Rh$Rs_comp <- "Rh"
  colnames(sub_Rh) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Ra <- subset(sdata, !is.na(Ra_root), select = c(StudyNumber, Ra_root, Tm))
  sub_Ra$Source <- "MGRsD_TA"
  sub_Ra$Rs_comp <- "Ra"
  colnames(sub_Ra) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Rh_Ra <- rbind(sub_Rh, sub_Ra)
  return(sub_Rh_Ra)
}


# colnames(DGRhRa)
RhRa_D <- function (sdata) {
  sub_Rh <- subset(sdata, !is.na(Rh_Norm), select = c(StudyNumber, Rh_Norm, Tsoil ))
  sub_Rh$Source <- "DGRsD_TS"
  sub_Rh$Rs_comp <- "Rh"
  colnames(sub_Rh) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Ra <- subset(sdata, !is.na(Ra_Norm), select = c(StudyNumber, Ra_Norm, Tsoil))
  sub_Ra$Source <- "DGRsD_TS"
  sub_Ra$Rs_comp <- "Ra"
  colnames(sub_Ra) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Rh_Ra <- rbind(sub_Rh, sub_Ra)
  return(sub_Rh_Ra)
}


# comparing Q10 of different season by Biome
Q10_Season_Biome <- function(sdata) {
  # estimate Tsoil by Tm
  sdata$TS <- ifelse( !is.na(sdata$TS), sdata$TS, 3.171523 + 0.801877*sdata$Tm )
  # get season information
  sdata$Season <- ifelse(sdata$Latitude>=0
                         , as.character(cut(sdata$Measure_Month, breaks = c(0, 2, 5, 8, 11, 13), labels = c("d.Winter", "a.Spring", "b.Summer", "c.Fall", "d.Winter")) )
                         , as.character(cut(sdata$Measure_Month, breaks = c(0, 2, 5, 8, 11, 13), labels = c("b.Summer", "c.Fall", "d.Winter", "a.Spring", "b.Summer")) )
                         )
  # OSH, SAV, SNO, URB merge together due to lack of enough data
  sdata$IGBP_FromPaper <- ifelse(sdata$IGBP_FromPaper=='OSH'|sdata$IGBP_FromPaper=='SAV'|sdata$IGBP_FromPaper=='SNO'|sdata$IGBP_FromPaper=='URB'
                                 , "Other", as.character(sdata$IGBP_FromPaper))
  
  g_plot <- ggplot( data = sdata, aes(TS, log(Rs_Norm)) ) + geom_point(col = "gray") + geom_smooth(method = "lm") +
    facet_grid(Season~IGBP_FromPaper)
    
  print(g_plot)
}

# calculate Q10 of Rs, Rh, and Ra
Q10_cal <- function (sdata) {
  # Using T_Del if T_soil not available
  sdata$Tsoil <- ifelse(is.na(sdata$Tsoil), 3.192786 + 0.771291 * sdata$Tm_Del, sdata$Tsoil)
  
  results <- data.frame()
  var_Study <- sort(unique(sdata$StudyNumber))

  for(i in 1 : length(var_Study) ) {
    
    MGRsD_ID <- sdata[sdata$StudyNumber == var_Study[i], ]
    var_SiteID <- sort(unique(MGRsD_ID$SiteID))
    
    for (j in 1:length(var_SiteID) ) {
      MGRsD_Site <- MGRsD_ID[MGRsD_ID$SiteID == var_SiteID[j], ]
      
      # if (nrow(MGRsD_Site) < 7) { next }
      
      Q10_Rs_Paper <- mean(MGRsD_Site$Q10_0_10_Rs, na.rm = T)
      Q10_Rh_Paper <- mean(MGRsD_Site$Q10_0_10_Rh, na.rm = T)
      Q10_Ra_Paper <- mean(MGRsD_Site$Q10_Ra, na.rm = T)
      obs <- subset(MGRsD_Site, !is.na(MGRsD_Site$RS_Norm) ) %>% nrow() 
      obs_Rh <- subset(MGRsD_Site, !is.na(MGRsD_Site$Rh_Norm) ) %>% nrow()
      obs_Ra <- subset(MGRsD_Site, !is.na(MGRsD_Site$Ra_Norm) ) %>% nrow()
      Biome <- MGRsD_Site$IGBP_FromPaper[1]
      Partition_Method <- MGRsD_Site$Partition_Method[1]
      
      
      # Calculate Q10_Rs and error handle
      m_nls <- try( nls(RS_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                        , data = MGRsD_Site, start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
      
      if(isTRUE(class(m_nls)=="try-error" | obs < 6 )) { 
          R <- NA
          p_R <- NA
          Q <- NA
          p_Q <- NA
          Mean_TS <- NA
          Q10_Rs <- NA
        } else { 
          sum_nls <- summary(m_nls)
          R <- sum_nls$coef[1, c(1)]
          p_R <- round(sum_nls$coef[1, c(4)], 4)
          Q <- sum_nls$coef[2, c(1)]
          p_Q <- round(sum_nls$coef[2, c(4)], 4)
          Mean_TS <- mean(MGRsD_ID$Tsoil)
          Q10_Rs <- ifelse(p_Q < 0.1, exp(Q*10), NA )
          Q10_Rs <- ifelse(Q10_Rs <= 10 & Q10_Rs >= 1, Q10_Rs, NA )
        } 
      
      # Calculate Q10_Rh and error handle
      Rh_nls <- try( nls(Rh_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                        , data = MGRsD_Site, start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
      if( isTRUE(class(Rh_nls)=="try-error" | obs_Rh < 6 ) ){
        Q_Rh <- NA
        p_Q_Rh <- NA
        Q10_Rh <- NA
      } else {
        sum_nls_Rh <- summary(Rh_nls)
        Q_Rh <- sum_nls_Rh$coef[2, 1]
        p_Q_Rh <- sum_nls_Rh$coef[2, 4] %>% round(4)
        Q10_Rh <- ifelse(p_Q_Rh < 0.1, exp(Q_Rh*10), NA)  
        Q10_Rh <- ifelse(Q10_Rh <= 10 & Q10_Rh >= 1, Q10_Rh, NA )
      }
      
      # Calculate Q10_Ra and error handle
      Ra_nls <- try( nls(Ra_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                         , data = MGRsD_Site, start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
      if( isTRUE(class(Ra_nls)=="try-error" | obs_Ra < 6) ){
        Q_Ra <- NA
        p_Q_Ra <- NA
        Q10_Ra <- NA
      } else {
        sum_nls_Ra <- summary(Ra_nls)
        Q_Ra <- sum_nls_Ra$coef[2, 1]
        p_Q_Ra <- sum_nls_Ra$coef[2, 4] %>% round(4)
        Q10_Ra <- ifelse(p_Q_Ra < 0.1, exp(Q_Ra*10), NA)
        Q10_Ra <- ifelse(Q10_Ra <= 10 & Q10_Ra >= 1, Q10_Ra, NA )
      }
      
        
      results <- rbind(results
                       ,data.frame(var_Study[i], var_SiteID[j], R, p_R, Q, p_Q, Mean_TS, Q10_Rs, obs
                                   , Q_Rh, p_Q_Rh, Q10_Rh, obs_Rh
                                   , Q_Ra, p_Q_Ra, Q10_Ra, obs_Ra
                                   , Biome, Partition_Method, Q10_Rs_Paper, Q10_Rh_Paper, Q10_Ra_Paper
                                   ) )
      
      print(paste0('*****', i,':' ,var_Study[i], '*****',j,":" ,var_SiteID[j]))
    }
  }
  
  # colnames(results) <- c( "StudyNumber", "SiteID", "R", "p_R", "Q", "p_Q", "Tm_Annual", "Q10_Rs", "obs", "Biome" )
  return (results)
}

#--------------------------------------------------------------------------------
# plot Rs, Rh, and Ra vs Tsoil by study ID and site
Rs_Rh_Ra_Tsoil <- function (sdata) {
  sdata$Tsoil <- ifelse(is.na(sdata$Tsoil), 3.192786 + 0.771291 * sdata$Tm_Del, sdata$Tsoil)
  
  var_Study <- sort(unique(sdata$StudyNumber))
  
  for(i in 1 : length(var_Study) ) {
    
    MGRsD_ID <- sdata[sdata$StudyNumber == var_Study[i], ]
    var_SiteID <- sort(unique(MGRsD_ID$SiteID))
    
    for (j in 1:length(var_SiteID) ) {
      MGRsD_Site <- MGRsD_ID[MGRsD_ID$SiteID == var_SiteID[j], ]
      
      # plot Rs_vs_Tsoil 
      plot(x=Tsoil, y=RS_Norm, data = sdata) 
      
      # plot Rh_vs_Tsoil
      points()
      
      # plot Ra_vs_Tsoil
      
      print(paste0('*****', i,':' ,var_Study[i], '*****',j,":" ,var_SiteID[j]))
    }
  }
}

#--------------------------------------------------------------------------------
# plot Rs Rh and Ra
Q10_comparison <- function (sdata) {
  
  # filt_data <- subset(sdata, sdata$p_Q < 0.1 & sdata$p_Q_Rh < 0.1 & sdata$p_Q_Ra < 0.1)
  # filt_data <- subset(filt_data, filt_data$Q10_Ra < 10)
  # filt_data <- subset(filt_data, filt_data$Q10_Rs > 1 & filt_data$Q10_Rh > 1 & filt_data$Q10_Ra > 1)
  # sdata <- subset(sdata, sdata$Q10_Rs_Paper < 10)
  
  title_Rs <- paste0("(a) Q10 of Rs, n=", nrow( subset(sdata, !is.na(sdata$Q10_Rs)))
                         , ", mean=", mean(sdata$Q10_Rs, na.rm = TRUE) %>% round(2) ) 
  title_Rs_Paper <- paste0("(b) Q10 of Rs(Paper), n=", nrow( subset(sdata, !is.na(sdata$Q10_Rs_Paper)) )
                                  , ", mean=", mean(sdata$Q10_Rs_Paper, na.rm = TRUE) %>% round(2) ) 
  
  
  plot_Rs <- ggplot(sdata, aes(Q10_Rs) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Rs) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
    # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  plot_Rs_Paper <- ggplot(sdata, aes(Q10_Rs_Paper) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Rs_Paper) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_Rs_Paper, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  
  # plot Q10 of Rh
  title_Rh <- paste0("(c) Q10 of Rh, n=", nrow( subset(sdata, !is.na(sdata$Q10_Rh)) )
                         , ", mean=", mean(sdata$Q10_Rh, na.rm = TRUE) %>% round(2) ) 
  title_Rh_Paper <- paste0("(d) Q10 of Rh(Paper), n=", nrow( subset(sdata, !is.na(sdata$Q10_Rh_Paper)) )
                         , ", mean=", mean(sdata$Q10_Rh_Paper, na.rm = TRUE) %>% round(2) ) 
  
  plot_Rh <- ggplot(sdata, aes(Q10_Rh) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Rh) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0,8) +
    geom_vline(xintercept = mean(sdata$Q10_Rh, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  
  plot_Rh_Paper <- ggplot(sdata, aes(Q10_Rh_Paper) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Rh_Paper) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0,8) +
    geom_vline(xintercept = mean(sdata$Q10_Rh_Paper, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  
  # plot Q10 of Ra
  title_Ra <- paste0("(e) Q10 of Ra, n=", nrow( subset(sdata, !is.na(sdata$Q10_Ra)) )
                         , ", mean=", mean(sdata$Q10_Ra, na.rm = TRUE) %>% round(2) ) 
  title_Ra_Paper <- paste0("(f) Q10 of Ra(Paper), n=", nrow( subset(sdata, !is.na(sdata$Q10_Ra_Paper)) )
                         , ", mean=", mean(sdata$Q10_Ra_Paper, na.rm = TRUE) %>% round(2) ) 
  
  plot_Ra <- ggplot(sdata, aes(Q10_Ra) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Ra) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0,8) +
    geom_vline(xintercept = mean(sdata$Q10_Ra, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  
  plot_Ra_Paper <- ggplot(sdata, aes(Q10_Ra_Paper) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_Ra_Paper) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0,8) +
    geom_vline(xintercept = mean(sdata$Q10_Ra_Paper, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  
  plot_grid(plot_Rs, plot_Rs_Paper, plot_Rh, plot_Rh_Paper, plot_Ra, plot_Ra_Paper, ncol = 2)
}

# agg <- aggregate(RS_Norm~StudyNumber + SiteID, FUN = length, data = DGRhRa )
# agg[order(agg$RS_Norm),]

check_outler <- function (sdata, study_number) {
  sub_data <- subset( sdata, sdata$StudyNumber == study_number )
  Rs_plot <- qplot( Tm_Del, RS_Norm, data = sub_data)
  Rh_plot <- qplot( Tm_Del, Rh_Norm, data = sub_data)
  Ra_plot <- qplot( Tm_Del, Ra_Norm, data = sub_data)
  plot_grid (Rs_plot, Rh_plot, Ra_plot, ncol = 3)
}


#--------------------------------------------------------------------------------
# plot Rh and Ra group by measure method

RhRa_partition <- function (sdata, sdata2) {
  sdata$Q10_Rh <- ifelse (!is.na(sdata$Q10_Rh_Paper), sdata$Q10_Rh_Paper, sdata$Q10_Rh)
  sdata$Q10_Ra <- ifelse (!is.na(sdata$Q10_Ra_Paper), sdata$Q10_Ra_Paper, sdata$Q10_Ra)
  
  Rh_P <- ggplot(sdata, aes(x=Partition_Method, y=Q10_Rh, group = Partition_Method)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  Ra_P <- ggplot(sdata, aes(x=Partition_Method, y=Q10_Ra, group = Partition_Method)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') + 
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  HC_P <- ggplot(sdata2, aes(x=Partition_Method, y=HC, group = Partition_Method)) + geom_violin() +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  RC_P <- ggplot(sdata2, aes(x=Partition_Method, y=RC, group = Partition_Method)) + geom_violin() +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
    
  plot_grid(Rh_P, Ra_P, HC_P, RC_P, ncol=1)
}


RhRa_biome <- function (sdata, sdata2) {
  sdata$Q10_Rh <- ifelse (!is.na(sdata$Q10_Rh_Paper), sdata$Q10_Rh_Paper, sdata$Q10_Rh)
  sdata$Q10_Ra <- ifelse (!is.na(sdata$Q10_Ra_Paper), sdata$Q10_Ra_Paper, sdata$Q10_Ra)
  
  Rh_P <- ggplot(sdata, aes(x=Biome, y=Q10_Rh, group = Biome)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  Ra_P <- ggplot(sdata, aes(x=Biome, y=Q10_Ra, group = Biome)) + geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') + 
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  HC_P <- ggplot(sdata2, aes(x=IGBP_FromPaper, y=HC, group = IGBP_FromPaper)) + geom_violin() +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  RC_P <- ggplot(sdata2, aes(x=IGBP_FromPaper, y=RC, group = IGBP_FromPaper)) + geom_violin() +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  plot_grid(Rh_P, Ra_P, HC_P, RC_P, ncol=1)
}


RhRa_month <- function (sdata2) {
  
  HC_P_IGBP <- ggplot(sdata2, aes(x=Measure_Month, y=HC, group = Measure_Month)) + geom_violin() +
    facet_grid(IGBP_FromPaper~.) +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  HC_P_Climate <- ggplot(sdata2, aes(x=Measure_Month, y=HC, group = Measure_Month)) + geom_violin() +
    facet_grid(Biome_M~.) +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  
  RC_P_IGBP <- ggplot(sdata2, aes(x=Measure_Month, y=RC, group = Measure_Month)) + geom_violin() +
    facet_grid(IGBP_FromPaper~.) +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  RC_P_Climate <- ggplot(sdata2, aes(x=Measure_Month, y=RC, group = Measure_Month)) + geom_violin() +
    facet_grid(Biome_M~.) +
    # geom_jitter(shape=16, position=position_jitter(0.2), col = 'gray') +
    stat_summary(fun.y=median, geom="point", size=2, color="red")
  
  print(HC_P_IGBP)
  print(HC_P_Climate)
  print(RC_P_IGBP)
  print(RC_P_Climate)
  # plot_grid(HC_P, RC_P, ncol=1)
}



#*****************************************************************************************************************
# functions for HGRsD
#*****************************************************************************************************************
# test daytime and night time Q10

HGRsD_DayNight_test <- function (sdata) {
  
  sdata <- sdata[order(sdata$MiddleClimate,sdata$StudyNumber,sdata$TimeLable),] 
  sub_day <- subset(sdata, sdata$DN == 'Day')
  sub_night <- subset(sdata, sdata$DN == 'Night')
  
  print( qplot(Tsoil, RS_Norm, data  = sdata) + facet_grid (.~DN) )
  
  # ------------------------------ calculate Q10 of day time
  day_nls <- nls(RS_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                    , data = sub_day, start = list(R = 0.5, Q = 0.035), trace = TRUE)
  
  Q_day <- summary(day_nls)$coef[2, c(1)]
  Q10_day <- exp(Q_day * 10)
  
  print( summary (day_nls) )
  print (paste0('Day time Q10 = ', Q10_day))
  
  # ------------------------------ calculate Q10 of night time
  night_nls <- nls(RS_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                 , data = sub_night, start = list(R = 0.5, Q = 0.035), trace = TRUE)
  
  Q_night <- summary(night_nls)$coef[2, c(1)]
  Q10_night <- exp(Q_night * 10)
  
  print( summary(night_nls) )
  print (paste0('Night time Q10 = ', Q10_night))
  
}

#********************************************************************
# by study_site
Q10_DayNight_byStudy <- function (sdata) {
  # sub_day <- subset(sdata, sdata$DN == 'Day')
  # sub_night <- subset(sdata, sdata$DN == 'Night')
  
  outputs <- data.frame()
  var_Study <- sort(unique(sdata$StudyNumber))
  
  for(i in 1 : length(var_Study) ) {
    
    HGRsD_ID <- sdata[sdata$StudyNumber == var_Study[i], ]
    var_SiteID <- sort(unique(HGRsD_ID$SiteID))
    
    for (j in 1:length(var_SiteID) ) {
      HGRsD_Site <- HGRsD_ID[HGRsD_ID$SiteID == var_SiteID[j], ]
      
      # meta information
      Season <- HGRsD_Site$Measure_Season[1]
      Month <- HGRsD_Site$Measure_Month[1]
      Biome <- HGRsD_Site$MiddleVegType[1]
      Climate <- HGRsD_Site$MiddleClimate[1]
      obs_day <- subset(HGRsD_Site, HGRsD_Site$DN == 'Day') %>% nrow()
      obs_night <- subset(HGRsD_Site, HGRsD_Site$DN == 'Night') %>% nrow()
      Mean_TS <- mean(HGRsD_ID$Tsoil, na.rm = T)
      
      # Calculate day time Q10 and error handle
      day_nls <- try( nls(RS_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                        , data = subset(HGRsD_Site, HGRsD_Site$DN == 'Day'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
      
      if(isTRUE(class(day_nls)=="try-error" | obs_day < 6 )) { 
        R <- NA
        p_R <- NA
        Q <- NA
        p_Q <- NA
        Q10_day <- NA
      } else { 
        sum_nls <- summary(day_nls)
        R <- sum_nls$coef[1, c(1)]
        p_R <- round(sum_nls$coef[1, c(4)], 4)
        Q <- sum_nls$coef[2, c(1)]
        p_Q <- round(sum_nls$coef[2, c(4)], 4)
        Q10_day <- ifelse(p_Q < 0.1, exp(Q*10), NA )
        Q10_day <- ifelse(Q10_day <= 10 & Q10_day >= 1, Q10_day, NA )
      } 
      
      # Calculate night time Q10 and error handle
      night_nls <- try( nls(RS_Norm ~ R  * exp(Q * Tsoil), nls.control(maxiter=100)
                            , data = subset(HGRsD_Site, HGRsD_Site$DN == 'Night'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
      if( isTRUE(class(night_nls)=="try-error" | obs_night < 6 ) ){
        Q_night <- NA
        p_Q_night <- NA
        Q10_night <- NA
      } else {
        sum_nls_night <- summary(night_nls)
        Q_night <- sum_nls_night$coef[2, 1]
        p_Q_night <- sum_nls_night$coef[2, 4] %>% round(4)
        Q10_night <- ifelse(p_Q_night < 0.1, exp(Q_night*10), NA)  
        Q10_night <- ifelse(Q10_night <= 10 & Q10_night >= 1, Q10_night, NA )
      }
      
      outputs <- rbind(outputs
                       ,data.frame(var_Study[i], var_SiteID[j]
                                   , Season, Month, Biome, Climate, Mean_TS
                                   , R, p_R, Q, p_Q, Q10_day, obs_day
                                   , Q_night, p_Q_night, Q10_night, obs_night
                       ) )
      
      print(paste0('*****', i,':' ,var_Study[i], '*****',j,":" ,var_SiteID[j]))
    }
  }
  
  # colnames(outputs) <- c( "StudyNumber", "SiteID", "R", "p_R", "Q", "p_Q", "Tm_Annual", "Q10_Rs", "obs", "Biome" )
  return (outputs)
}


# Q10 of day night comparison
# plot Rs Rh and Ra
Q10_DN_Comp <- function (sdata) {
  
  title_day <- paste0("(a) Day time Q10, n=", nrow( subset(sdata, !is.na(sdata$Q10_day)) )
                     , ", mean=", mean(sdata$Q10_day, na.rm = TRUE) %>% round(2) ) 
  
  title_night <- paste0("(b) Night time Q10, n=", nrow( subset(sdata, !is.na(sdata$Q10_night)) )
                      , ", mean=", mean(sdata$Q10_night, na.rm = TRUE) %>% round(2) ) 
  
  # plot day time Q10
  plot_Q10_day <- ggplot(sdata, aes(Q10_day) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_day) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_day, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  # plot night time Q10
  plot_Q10_night <- ggplot(sdata, aes(Q10_night) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_night) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_night, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  plot_grid(plot_Q10_day, plot_Q10_night, ncol = 2)
}

# colnames(Q10_DN_site)

Q10_DN_Biome <- function (sdata) {
  
  day_sub <- subset (sdata, select = c("var_Study.i.",  "var_SiteID.j.", "Season", "Month", "Biome"        
                                       ,"Climate", "Mean_TS"
                                       , "Q","p_Q", "Q10_day", "obs_day") )
  
  colnames (day_sub) <- c("StudyID",  "SiteID", "Season", "Month", "Biome"        
                          ,"Climate", "Mean_TS", "Q", "p_Q", "Q10", "obs")
  day_sub$DN <- 'Day'
  
  night_sub <- subset(sdata, select = c("var_Study.i.",  "var_SiteID.j.", "Season", "Month", "Biome"        
                                        ,"Climate", "Mean_TS", "Q_night", "p_Q_night", "Q10_night", "obs_night")  )
  
  colnames (night_sub) <- c("StudyID",  "SiteID", "Season", "Month", "Biome"        
                          ,"Climate", "Mean_TS", "Q", "p_Q", "Q10", "obs")
  night_sub$DN <- 'Night'
  
  comb_data <- rbind(day_sub, night_sub)
  
  # plot Q10 of day and night facet by Biome
  plot_Q10_byBiome <- ggplot(comb_data, aes(Q10) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    facet_grid(Biome ~ DN)
  
  print(plot_Q10_byBiome)
}



