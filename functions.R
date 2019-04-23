
# functions

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
  sub_Rh <- subset(sdata, !is.na(Rh), select = c(StudyNumber, Rh, Tsoil.C. ))
  sub_Rh$Source <- "DGRsD_TS"
  sub_Rh$Rs_comp <- "Rh"
  colnames(sub_Rh) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Ra <- subset(sdata, !is.na(Ra_root), select = c(StudyNumber, Ra_root, Tsoil.C.))
  sub_Ra$Source <- "DGRsD_TS"
  sub_Ra$Rs_comp <- "Ra"
  colnames(sub_Ra) <- c("StudyNumber", "Rs", "Tm", "Source", "Rs_comp")
  
  sub_Rh_Ra <- rbind(sub_Rh, sub_Ra)
  return(sub_Rh_Ra)
}



