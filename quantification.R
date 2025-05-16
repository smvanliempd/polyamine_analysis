d.qnt <- d.all[Label_enrichment == 0
              ][Sample.Class == "Curve",
                 .(Assay.Date,
                   Run.Index,
                   Run = paste0("Run ", Run.Index),
                   Analyte = Metabolite,
                   Injection.Replicate,
                   Conc = Curve.Concentrations,
                   Sig = Signal_MFC)
              ][ , LLOQ :=12*mean(Sig[Conc == 0]), by = .(Assay.Date, Analyte)
              ][Conc > 0 & Sig > LLOQ]

d.qnt.pars <- d.qnt[, {
  d <- .SD
  d_low <- d[Conc <= 0.5 ]
  m_low <- lm(log(d_low$Sig) ~ log(d_low$Conc))
  summ_low <- summary(m_low)
  coefs_low <- summ_low$coefficients[1:2,1]
  r_adj_low <- summ_low$adj.r.squared
  
  d_high <- d[Conc >= 0.5 ]
  m_high <- lm(log(d_high$Sig) ~ log(d_high$Conc))
  summ_high <- summary(m_high)
  coefs_high <- summ_high$coefficients[1:2,1]
  r_adj_high <- summ_high$adj.r.squared
  
  lloq <- LLOQ[1]
  
  list(a_low  = coefs_low[1],
       b_low  = coefs_low[2], 
       r_low  = r_adj_low,
       x_low  = min(d$Conc),
       y_low  = max(d$Sig),
       
       a_high = coefs_high[1],
       b_high = coefs_high[2], 
       r_high  = r_adj_high,
       x_high  = max(d$Conc),
       y_high  = min(d$Sig),
       
       lloq = lloq,
       signal_cut = mean(Sig[Conc == 0.5])
       
       )
}, by = .(Assay.Date,Run, Analyte) ]

d.qnt <- d.qnt[d.qnt.pars[,.(a_low,b_low,b_high,a_high,signal_cut,Assay.Date,Run,Analyte)], on = .(Assay.Date, Analyte) ]
d.qnt[ , S_mod_low := exp(b_low*log(Conc) + a_low)
    ][ , S_mod_high := exp(b_high*log(Conc) + a_high)
    ][ , Conc_calc :=  ifelse(Sig < signal_cut, 
                              calc_conc(Sig,a_low,b_low), 
                              calc_conc(Sig,a_high,b_high))
    ][ , dev_from_theor := (Conc_calc/Conc) - 1 
    ][ , flag := ifelse(abs(dev_from_theor) > 0.15, "*","")]

# get and plot mean deviations ± sem from nominal concentrations
d.dev <- d.qnt[ , .(dev_mn =  round(100*mean(dev_from_theor),1),
                    dev_sem =  round(100*sd(dev_from_theor,1)/sqrt(max(d.qnt$Run.Index)),1)), by = .(Conc, Analyte)
             ][ , c("v_min","v_max","tab_val") := { .(v_min = dev_mn - dev_sem,
                                                      v_max = dev_mn + dev_sem,
                                                      tab_val = paste(dev_mn,"±",dev_sem)) } ]
tab.devs <- dcast(d.dev,Conc~Analyte, value.var = "tab_val")
p.dev <- ggplot(d.dev, aes(x = as.factor(Conc), y = dev_mn, col = Analyte))+
  geom_hline(yintercept = c(-15,0,15), lty = 2, col = "grey50")+
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = v_min,ymax = v_max),
                 position = position_dodge(width = 0.5))+
  scale_color_manual(name = "Analyte",values = c("#00baff","#ff4500","#ffc500"))+
  scale_y_continuous(breaks = seq(-50,50,5),
                     minor_breaks = seq(-50,50,1),guide = guide_axis(minor = TRUE))+
  labs(x = "Concentration (µM)", y = "%Deviation from nominal\n(mean ± sem)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        # legend.position = "bottom",
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.dev

p.qnt01 <- ggplot(d.qnt,
       aes(x = Conc))+
  geom_vline(xintercept = 0.5, lty = 2)+
  geom_hline(data = d.qnt.pars, aes(yintercept = lloq), lty =2 )+
  geom_text(data = d.qnt.pars, aes(x = 0.03, 
                                   y = y_low, 
                                   label = paste("R","[adj]","^2","==", 
                                                 format(round(r_low,3), nsmall = 0))), 
            parse = TRUE, col = "blue",size = 2)+
  geom_text(data = d.qnt.pars, aes(x = 3.5, 
                                   y = 6*y_high, 
                                   label = paste("R","[adj]","^2","==", 
                                                 format(round(r_high,3), nsmall = 0))), 
            parse = TRUE, col = "orangered",size = 2)+
  geom_point(aes(y = Sig, 
                 shape = factor(Injection.Replicate)),
             alpha = 0.6, size = 1)+
  geom_line(aes(y = S_mod_low), col = "blue", lty =2) +
  geom_line(aes(y = S_mod_high), col = "orangered", lty =2) +
  scale_y_continuous(transform = "log10",guide = guide_axis(minor = TRUE))+
  scale_x_continuous(transform = "log10",guide = guide_axis(minor = TRUE))+
  scale_shape_manual(name = "Curve\nnumber", values = c(16,17))+
  facet_grid( Analyte~Run , scales = "free_y")+
  labs(x = "Concentration (µM)", y = "Signal")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        legend.position = "bottom",
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.qnt01

run.sel <- "Run 4"
p.qnt02 <- ggplot(d.qnt[ Run == run.sel],
                  aes(x = Conc))+
  geom_vline(xintercept = 0.5, lty = 2)+
  geom_hline(data = d.qnt.pars[Run== run.sel], aes(yintercept = lloq), lty =2 )+
  geom_text(data = d.qnt.pars[Run== run.sel], aes(x = 0.03, 
                                   y = y_low, 
                                   label = paste("R","[adj]","^2","==", 
                                                 format(round(r_low,3), nsmall = 0))), 
            parse = TRUE, col = "blue",size = 2)+
  geom_text(data = d.qnt.pars[Run== run.sel], aes(x = 3.5, 
                                   y = 6*y_high, 
                                   label = paste("R","[adj]","^2","==", 
                                                 format(round(r_high,3), nsmall = 0))), 
            parse = TRUE, col = "orangered",size = 2)+
  geom_point(aes(y = Sig, 
                 shape = factor(Injection.Replicate)),
             alpha = 0.6, size = 1)+
  geom_line(aes(y = S_mod_low), col = "blue", lty =2) +
  geom_line(aes(y = S_mod_high), col = "orangered", lty =2) +
  scale_y_continuous(transform = "log10",guide = guide_axis(minor = TRUE))+
  scale_x_continuous(transform = "log10",guide = guide_axis(minor = TRUE))+
  scale_shape_manual(name = "Curve\nnumber", values = c(16,17))+
  facet_wrap(.~ Analyte , scales = "free_y")+
  labs(x = "Concentration (µM)", y = "Signal")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title  = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.qnt02

