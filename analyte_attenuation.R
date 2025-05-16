# get data
d.isp <-d.all[(Experiment %in% c("Stability","Use case") & Issup == "x"),
              .(Analyte = Metabolite,
                Assay.Date,
                Run.Index, 
                Injection.ID, 
                Signal,
                Signal_MFC_QC,
                Conc_raw) 
            ][ , Group := ifelse(grepl("Curve", Injection.ID), "SOL", ifelse(grepl("QC", Injection.ID),"NO", "SPIKE"))
            ][ , grp := ifelse(grepl("Curve", Injection.ID), 3, ifelse(grepl("QC", Injection.ID),1, 2))
            ][ , spl := str_extract(Injection.ID, "\\d(?=[^\\d]*$)")
            ][ , run := as.integer(factor(Run.Index)  )]
setkeyv(d.isp, c("run","grp","spl"))
d.isp <- get_conc(dt = d.isp, signal_type = "Signal")

# fit
mod_isupp_04B <- stan_model("isupp_04B.stan")
fits.issup <- sapply(mets, function(m) {
  dd <- d.isp[ Analyte == m & !(run %in% c(1,3))
            ][ , run := as.integer(factor(run))]
  dstan <- list(N = nrow(dd),
                N_grp = max(dd$grp),
                N_run = max(dd$run),
                grp = dd$grp,
                run =  dd$run ,
                Y = dd$Signal)
  fit <- sampling(mod_isupp_04B, data = dstan,
                  chains =1,  cores=4,
                  warmup = 3000, iter = 6000, 
                  control = list(adapt_delta = 0.99, max_treedepth = 12),
                  seed = seed)
  return(fit)
}, USE.NAMES = TRUE)

# extract samples 
smpl.isp <-  sapply(mets, function(m){
  mod <- fits.issup[[m]]
  ss <- extract(mod,pars = c("S_1000")) |> 
    do.call(cbind, args= _) |>
    data.table() |>
    melt(data = _,measure.vars=1)
  ss <- ss[ , .(S = quantile(value,qtl_vec),
                qntl = paste0("Q",round(100*qtl_vec,0)),
                Analyte = m), by = variable]
  ss <- dcast(ss, ...~qntl, value.var = "S")
},simplify = FALSE)
smpl.isp <- do.call(rbind, smpl.isp)

#  plot
p.isp <- ggplot(smpl.isp, aes(x = Analyte, y = Q50, col = variable ) ) +
  geom_point(position = position_dodge(width = .5), size = 1) +
  geom_linerange(aes(ymin = Q10, ymax = Q90),position = position_dodge(width = .5), linewidth = 0.25)+
  geom_linerange(aes(ymin = Q25, ymax = Q75),position = position_dodge(width = .5), linewidth = 1)+
  scale_color_manual(name = "Spike level ",values = c("#ff4500","#00baff","#ffc500"), labels = c("1 µM"))+
  labs(y = "Signal attenuation (median, QI50/80)")+
  scale_y_continuous(breaks = seq(0,1,0.2),
                     minor_breaks = seq(0,1,0.05),
                     guide = guide_axis(minor = TRUE),
                     limits = c(0,NA))+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.8,0.15),
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
p.isp

# ggplot(d.isp, aes( x = Group, y = Conc_raw_alt, col = factor(run) ))+
#   geom_point(position = position_dodge(width = 0.5))+
#   facet_grid(Analyte~., scales = "free_y")+
#   ylim(0,NA)+
#   theme_bw() 

