# get/clean data
d.rec <-  d.all[Experiment == "Recovery/Ion suppresion" & Sample.Type != "SOL" ,
            .(Analyte = Metabolite,
              Sample.Type,
              Concentration,
              spl     = Sample.Index, 
              grp     = Group.Index,
              Conc_raw)]
setkeyv(d.rec,c("grp","spl"))

# fit
mod0_recov_06A <- stan_model("recov_06A.stan")
fits.recov <- sapply(mets, function(m) {
  dd <- d.rec[Analyte == m]
  dstan <- list(N = nrow(dd),
                N_grp = max(dd$grp),
                N_spl = max(dd$spl),
                spl = dd$spl,
                grp = dd$grp,
                mu_grp = dd[ , mean(Conc_raw), by = grp]$V1,
                Y = dd$Conc_raw)
  fit <- sampling(mod0_recov_06A, data = dstan,
                  chains =1,  cores=4 ,
                  warmup = 3000, iter = 6000, 
                  control = list(adapt_delta = 0.99, max_treedepth = 12),
                  seed = seed)
  return(fit)
}, USE.NAMES = TRUE)

# extract samples
pars_vec <- c("R_100","R_1000")
smpl.rec <- sapply(mets, function(m){
  mod <- fits.recov[[m]]
  ss <- extract(mod,pars = pars_vec) |>
    do.call(cbind, args= _) |>
    data.table() |>
    melt(data = _,measure.vars=1:length(pars_vec))
  ss <- ss[ , .(R = 100*quantile(value,qtl_vec),
                qntl = paste0("Q",round(100*qtl_vec,0)),
                Analyte = m), by = variable]
  ss <- dcast(ss, ...~qntl, value.var = "R")
},simplify = FALSE)
smpl.rec <- do.call(rbind, smpl.rec)

# plot
p.rec <- ggplot(smpl.rec, aes(x = Analyte, y = Q50, col = variable ) ) +
  geom_hline(yintercept = c(80,100), lty = 2,col = "grey50")+
  geom_point(position = position_dodge(width = .5), size = 1) +
  geom_linerange(aes(ymin = Q10, ymax = Q90),position = position_dodge(width = .5), linewidth = 0.25)+
  geom_linerange(aes(ymin = Q25, ymax = Q75),position = position_dodge(width = .5), linewidth = 0.75)+
  scale_color_manual(name = "Spike level ",
                     values = c("#00baff","#ff4500","#ffc500"),
                     labels = c("100 nM", "1 µM"))+
  labs(y = "%Recovery (median, QI50/80)")+
  scale_y_continuous(breaks = 100*seq(0,1,0.2),
                     minor_breaks = 100*seq(0,1,0.05),
                     guide = guide_axis(minor = TRUE),
                     limits = c(0,NA))+
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
        plot.background = element_rect(colour = "white"),
        plot.margin = unit(c(0,0,0,30),"pt" # margin for patchwork
                           ))
p.rec

