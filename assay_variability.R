# get/clean data
d.var <- d.all[Sample.Group %in% c("Cells","QC")[1] & Experiment == "Stability"
             ][,.(Analyte = Metabolite ,
                  Assay.Date,
                  Run.Index,
                  Technical.Replicate,
                  Sample.ID,
                  Sample.Class,
                  Signal, 
                  Signal_MFC_QC,
                  Conc_raw)
             ][ , spl := Technical.Replicate + 6*(Run.Index-1)]
setkeyv(d.var, c("Run.Index","spl"))
d.var <- get_conc(dt=d.var,signal_type = "Signal_MFC_QC") 
d.var[Analyte == "Putrescine" & Run.Index == 1,c("Conc_raw_alt","Conc_raw") := NA]

# plot raw data
p.var01 <- ggplot(d.var, 
       aes(x = Run.Index, y = Conc_raw_alt,  grp = factor(spl)))+
  geom_hline(yintercept = 0,col = "black")+
  geom_hline(data = d.var[,.(mu = mean(Conc_raw,na.rm = TRUE)),by = Analyte],
             aes(yintercept = mu)
             , lty =2,col = "grey50" )+
  geom_point(position = position_dodge(width = 0.4), 
             size = 1,alpha = 0.4, show.legend = FALSE)+
  scale_x_continuous(breaks = 1:4, labels = paste0("day ",1:4))+
  facet_grid(Analyte~.,scales = "free")+
  labs( y = "Concentration in vial (µM)")+
  scale_y_continuous(breaks = seq(0,2,0.2),
                     minor_breaks = seq(0,2,0.05),guide = guide_axis(minor = TRUE), limits = c(0,NA))+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.var01
ggsave("Fig06_Stability_raw.png",plot = p.var01,path =paste0(getwd(),"/graphs_tables"),width = 5, height = 5,dpi = 600 )

# fit
var_mod04 <- stan_model("var_meanscld_04.stan") #non-centered version of 3, samples better than 3, used for analysis
fits.var <- sapply(mets, function(m) {
  if (m == "Putrescine") {
    # to delete first run for putrescine
    d <- d.var[Analyte == m & Run.Index != 1
           ][ ,Run.Index := Run.Index - 1
           ][ ,Spl := Technical.Replicate + 6*(Run.Index-1)]
  } else {
    d <- d.var[Analyte == m ]
  }
  dstan <- list(N= nrow(d),
                Y = d$Conc_raw,
                N_run = max(d$Run.Index),
                N_spl = max(d$spl),
                run = d$Run.Index,
                spl = d$spl)
  fit <- sampling(var_mod04, data = dstan,
                  chains =1,  cores = 4,
                  warmup = 3000, iter = 6000, init = 0,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  seed = seed)
  return(fit)
}, simplify = FALSE, USE.NAMES = TRUE)

# extract samples
par_names_var <- c("sigma", "sigma_spl", "sigma_run")
smpl.var <- sapply(mets, function(m){
  s <- extract(fits.var[[m]], pars = par_names_var)
  s <- do.call(cbind,s) |>
    data.table() |>
    melt(data = _, measure.vars = 1:length(par_names_var)) 
  ss <- s[,.(value = 100*quantile(value,qtl_vec), 
             qtl = paste0("Q",round(100*qtl_vec,0))),by = variable]
  ss$Analyte <- m
  ss <- dcast(ss,...~qtl, value.var = "value")
  return(ss)
}, simplify = FALSE)
smpl.var <- do.call(rbind, smpl.var)
smpl.var[ , metric := ifelse(variable == "sigma","measurement", ifelse(variable== "sigma_spl", "intra-day", "inter-day"))]
smpl.var[ , metric :=  factor(metric, c("inter-day","intra-day","measurement"))]

# plot model output
p.var02 <- ggplot(smpl.var, aes(x = Analyte, y = Q50)) +
  geom_hline(yintercept = c(5, 15), lty = 2, col = "grey50")+
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Q25, ymax = Q75), linewidth = 0.75)+
  geom_linerange(aes(ymin = Q10, ymax = Q90), linewidth = 0.25)+
  facet_wrap(.~metric,scales = "free")+
  labs(y = "%CV (median, QI50/80)")+
  scale_y_continuous(breaks = seq(0,100,5),
                     minor_breaks = seq(0,100,1),guide = guide_axis(minor = TRUE), limits = c(0,31))+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(linewidth  = 0.5,inherit.blank = FALSE),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.var02

# calculate minimal detectable effects
par_names_mde <- c("sigma_pooled_run","sigma_pooled_spl")
calculate_mde_alpha <- function(sigma_pooled, n, alpha = 0.05) {
  df <- 2 * n - 2                   # Degrees of freedom for a two-sample t-test
  t_alpha <- qt(1 - alpha / 2, df)  # Two-tailed critical t-value
  SE <- sigma_pooled * sqrt(2 / n)  # Standard error of the difference
  MDE <- t_alpha * SE               # Minimal Detectable Effect
  
  return(MDE)
}
n_mde <- 6
smpl.mde <- sapply(mets, function(m){
  s <- extract(fits.var[[m]], pars = par_names_mde)
  s <- do.call(cbind,s) |>
    data.table() |>
    melt(data = _, measure.vars = 1:length(par_names_mde)) 
  ss <- s[,.(value = calculate_mde_alpha(sigma_pooled = 100*quantile(value,qtl_vec),
                                         n = n_mde), 
             qtl = paste0("Q",round(100*qtl_vec,0))),by = variable]
  ss$Analyte <- m
  ss <- dcast(ss,...~qtl, value.var = "value")
  return(ss)
}, simplify = FALSE)
smpl.mde <- do.call(rbind, smpl.mde)
smpl.mde[ , metric := ifelse(variable == "sigma_pooled_run",
                                            paste0("mimal effect\nbetween run\n(n = ",n_mde,")"), 
                                            paste0("mimal effect\nwithin run\n(n = ",n_mde,")"))]


p.var03 <- ggplot(smpl.mde, aes(x = Analyte, y = Q50)) +
  geom_hline(yintercept = c(10, 35), lty = 2, col = "grey50")+
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Q25, ymax = Q75), linewidth = 0.75)+
  geom_linerange(aes(ymin = Q10, ymax = Q90), linewidth = 0.25)+
  facet_wrap(.~metric,scales = "free")+
  # labs(y = "sigma (QI50/80, mu_data = 1)")+
  labs(y = expression(paste("%change from H"[0], " (median, QI50/80)")))+
  scale_y_continuous(breaks = seq(0,100,5),
                     minor_breaks = seq(0,100,1),guide = guide_axis(minor = TRUE), limits = c(0,41))+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(linewidth  = 0.5,inherit.blank = FALSE),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.var03

# 

# ggplot(d.var, 
#        aes(x = Run.Index, grp = factor(spl)))+
#   # geom_hline(data = d.var[,.(mu = mean(Conc_raw_alt,na.rm = TRUE)),by = Analyte],
#   #            aes(yintercept = mu)
#   #            , lty =2,col = "grey50" )+
#   geom_point(aes(y = Conc_raw),position = position_dodge(width = 0.4), 
#              size = 1,alpha = 0.4, show.legend = FALSE)+
#   geom_point(aes(y = Conc_raw_alt),position = position_dodge(width = 0.4), 
#              size = 1,alpha = 0.4, show.legend = FALSE,col ="red")+
#   scale_x_continuous(breaks = 1:4, labels = paste0("day ",1:4))+
#   facet_grid(Analyte~.,scales = "free")+
#   labs( y = "Concentration in vial (µM)")+
#   scale_y_continuous(breaks = seq(0,2,0.2),
#                      minor_breaks = seq(0,2,0.05),guide = guide_axis(minor = TRUE), limits = c(0,NA))+
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
#         # axis.title.y = element_text(size = 8),
#         axis.line = element_line(colour = "grey50"),
#         axis.ticks = element_line(colour = "black"),
#         axis.ticks.length = unit(0.1,"cm"),
#         axis.minor.ticks.length = rel(0.5),
#         plot.background = element_rect(colour = "white"))
