# get data
media_types <- c("DMEM","Dialyzed","HPLM")
d.use.med <- d.all[Experiment == "Use case" & Sample.Class == 'Sample' & Sample.Group != "Spike" & Run.Index == 5,
               .(Analyte = Metabolite ,
                 Assay.Date,
                 grp = Group.Index,
                 spl = Sample.Index,
                 Injection.Replicate,
                 Sample.Group = factor(Sample.Group, media_types),
                 Signal,
                 Signal_MFC_QC,
                 Conc_raw)]
d.use.med <- get_conc(dt=d.use.med, signal_type = "Signal_MFC_QC")

# fit
mod_groups_01 <- stan_model("group_analysis_01.stan")
fits.use.med <- sapply(mets, function(m) {
  dd <- d.use.med[Analyte == m]
  dstan <- list(Y = dd$Conc_raw,
                spl = dd$spl,
                grp = dd$grp,
                N = nrow(dd),
                N_spl = max(dd$spl),
                N_grp = max(dd$grp)
  )
  fit <- sampling(mod_groups_01 , data = dstan,
                  chains =1,  cores=4 ,
                  warmup = 3000,
                  iter = 6000,
                  control = list(adapt_delta = 0.99, max_treedepth = 12),
                  seed = seed)
  return(fit)
}, USE.NAMES = TRUE)

# extract samples
smpl.use.med <- sapply(mets, function(m){
  mod <- fits.use.med[[m]]
  ss <- extract(mod,pars = c("s_grp")) |>
    do.call(cbind, args= _) |>
    data.table()
  setnames(ss, c("V1","V2","V3"), media_types)
  ss <-   melt(data = ss, measure.vars=1:3)
  ss <- ss[ , .(R = quantile(value,qtl_vec),
                qntl = paste0("Q",round(100*qtl_vec,0)),
                Analyte = m), by = variable]
  ss <- dcast(ss, ...~qntl, value.var = "R")
},simplify = FALSE)
smpl.use.med <- do.call(rbind, smpl.use.med)
smpl.use.med[ , variable := factor(variable, media_types)
       ][ , grp := as.integer(variable)]
d.use.med.agr <- d.use.med[ ,.(conc_mean = mean(Conc_raw)), by = .(Analyte,Sample.Group,spl,grp) ]

# plot
p.use.med <- ggplot(smpl.use.med, aes(x = grp - 0.2, y = Q50)) +
  geom_point(size =  1)+
  geom_linerange(aes(ymin = Q25, ymax = Q75), linewidth = 0.75)+
  geom_linerange(aes(ymin = Q10, ymax = Q90), linewidth = 0.25)+
  geom_point(data = d.use.med.agr, aes(x = grp +0.2,y = conc_mean),shape = 21, fill = "#FF000070",size =1)+
  facet_wrap(Analyte~.,scales = "free")+
  labs(x = "Media types", y = "Concentration in vial (µM)")+
  scale_x_continuous(breaks = 1:3, labels = media_types, limits = c(0.5,3.5))+
  scale_y_continuous(breaks = seq(0,3,0.2),
                     minor_breaks = seq(0,3,0.05),guide = guide_axis(minor = TRUE), limits = c(0,NA))+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.line.y = element_line(linewidth  = 0.5,inherit.blank = FALSE),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.use.med
ggsave("Fig09_culture_media_analysis.png",plot = p.use.med,path =paste0(getwd(),"/graphs_tables"),width = 5, height = 5,dpi = 600 )

