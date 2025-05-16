# get/clean data
d.ext <- d.all[Experiment == "Extraction" , 
               .(Analyte  = Metabolite,
                 Smpl = Sample.ID,
                 Injc = Injection.ID,
                 Repl = Injection.Replicate,
                 Ext  = Extraction.Method,
                 spl = Sample.Index,
                 grp = Group.Index,
                 run  = Run.Index,
                 Sig  = Signal_MFC_QC)
][ , Sig_scl := Sig/mean(Sig), by = .(Analyte)]
setkeyv(d.ext, c("run","grp","spl"))

# fit
mod_extract_05 <- stan_model("extract_05.stan")
fits.ext <- sapply(mets, function(m) {
  dd <- d.ext[Analyte == m]
  dstan <- list(Y_sm = dd$Sig_scl,
                spl = dd$spl,
                grp = dd$grp,
                run = dd$run,
                N = nrow(dd),
                N_run = max(dd$run),
                N_spl = max(dd$spl),
                N_grp = max(dd$grp)
                )
  fit <- sampling(mod_extract_05 , data = dstan,
                  chains =1,  cores=4 ,
                  warmup = 3000,
                  iter = 6000,
                  control = list(adapt_delta = 0.99, max_treedepth = 12),
                  seed = seed)
  return(fit)
}, USE.NAMES = TRUE)

# extract samples
smpl.ext <- sapply(mets, function(m){
  mod <- fits.ext[[m]]
  ss <- extract(mod,pars = c("a_grp")) |>
    do.call(cbind, args= _) |>
    data.table()
  setnames(ss, c("V1","V2","V3"), c("M80","M50","BP"))
  ss <- melt(data = ss, measure.vars=1:3)
  ss <- ss[ , .(R = quantile(value,qtl_vec),
                qntl = paste0("q",round(100*qtl_vec,0)),
                Analyte = m), by = variable]
  ss <- dcast(ss, ...~qntl, value.var = "R")
},simplify = FALSE)
smpl.ext <- do.call(rbind, smpl.ext)

# plot
p.ext <- ggplot(smpl.ext, aes(x = Analyte, y = q50, col = variable ) ) +
  geom_hline(yintercept = 1, lty =2,col = "grey50")+
  geom_point(position = position_dodge(width = .5), size = 1) +
  geom_linerange(aes(ymin = q10, ymax = q90),position = position_dodge(width = .5), linewidth = 0.25)+
  geom_linerange(aes(ymin = q25, ymax = q75),position = position_dodge(width = .5), linewidth = 0.75)+
  scale_color_manual(name = "Extraction\n method",values = c("#00baff","#ff4500","#ffc500"))+
  labs(y = "Scaled signal (median in QI50/80)")+
  scale_y_continuous(minor_breaks = seq(0,2,0.1),guide = guide_axis(minor = TRUE), limits = c(0,NA))+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.2,0.15),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1,"inch"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.ext

