lmer.diag.mplot = function(lmer.mod){

  #set up data frame with necessary vectors for plotting
  mod.dat = as.data.frame(as.numeric(fitted(lmer.mod)))
  names(mod.dat) = "fitted"
  mod.dat$resid = as.numeric(resid(lmer.mod)) #stanardized residuals
  mod.dat$qq.line.y = quantile(mod.dat$resid[!is.na(mod.dat$resid)], c(0.25, 0.75))
  mod.dat$qq.line.x = qnorm(c(0.25, 0.75))
  mod.dat$hatvalues = hatvalues(lmer.mod)
  
  # I. residuals vs. fitted plot
  resid.vs.fitted.plot = 
    ggplot(mod.dat, aes(x=fitted, y=resid)) +
    geom_point(shape=1) + 
    geom_smooth(color="red", size=0.5) +
    labs(title = "A. Model residuals vs. fitted values", 
         x = "fitted values", 
         y = "standardized residuals") +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  #II. Histogram for checking normality of residuals
  resid.hist = 
    ggplot(mod.dat, aes(x=resid)) + 
    geom_histogram(color="black", fill="white", bins=30) +
    labs(title = "B. Histogram of model residuals", 
         x = "standardized residuals") +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  # III. QQ-plot
  #calculate qqline for plotting
  slope = diff(mod.dat$qq.line.y)/diff(mod.dat$qq.line.x)
  int = mod.dat$qq.line.y[1L] - slope * mod.dat$qq.line.x[1L]
  #QQ-plot
  qq.plot = 
    ggplot(mod.dat, aes(sample = resid)) + 
    stat_qq(shape=1) + 
    geom_abline(slope = slope, intercept = int, size=0.5, color="red") +
    labs(title = "C. Quantile-quantile plot", 
         x = "theoretical quantiles", 
         y = "sample quantiles") +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
  
  #IV. Checking for leverage points
  lev.plot = 
    ggplot(mod.dat, aes(x=hatvalues , y=resid)) + 
    geom_point(shape=1) +
    geom_smooth(color="red", size=0.5) +
    labs(title = "D. Leverage plot", x = "leverage", y = "standardized residuals") +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))

  #create multiplot
  grid.arrange(resid.vs.fitted.plot, resid.hist, qq.plot, lev.plot, ncol = 2)
  
}



