# get data
inhib_types <- c("Veh","SAM486A")
d.use.inh <- d.all[Experiment == "Use case" & Sample.Class == 'Sample' & Sample.Group != "Spike" & Run.Index == 6,
                   .(Analyte = Metabolite ,
                     Lab =Label_enrichment,
                     Assay.Date,
                     grp = Group.Index,
                     spl = Sample.Index,
                     rep = Injection.Replicate,
                     Treatment = factor(Treatment, inhib_types),
                     Signal,
                     Signal_MFC_QC,
                     Conc_raw)]
d.use.inh <- get_conc(dt=d.use.inh, signal_type = "Signal_MFC_QC")
dt <- d.use.inh[!is.na(Conc_raw_alt)]
dtw <- dcast(d.use.inh,Treatment+Analyte+spl+grp+rep~Lab,value.var = "Conc_raw_alt")
dtw[,tot := sum(.SD,na.rm = TRUE), .SDcols = 6:8, by = 1:nrow(dtw)]
dtl <- melt(dtw,measure.vars = 6:9,value.name = "Conc_raw_alt",variable.name = "Lab",na.rm = TRUE)
dtl[,Analyte_lab := paste0(Analyte,Lab)]

mod_groups_inhib_02 <- stan_model("group_analysis_inhib_02.stan") 
mets_lab <- c("Putrescine0","Spermidine0","Spermidine3","Spermine0","Spermine3","Spermine6",
              "Putrescinetot","Spermidinetot","Sperminetot")
fits.use.inh <- sapply(mets_lab, function(m) {
    dd <- dtl[Analyte_lab == m ]
    dstan <- list(Y = dd$Conc_raw_alt,
                  spl = dd$spl,
                  grp = dd$grp,
                  N = nrow(dd),
                  N_spl = max(dd$spl),
                  N_grp = max(dd$grp)
    )
    fit <- sampling(mod_groups_inhib_02 , data = dstan,
                    chains =1,  cores=4 ,
                    warmup = 3000,
                    iter = 6000,
                    control = list(adapt_delta = 0.99, max_treedepth = 12),
                    seed = seed)
    return(fit)
}, USE.NAMES = TRUE)

# extract samples
smpl.use.inh <- sapply(mets_lab, function(m){
  mod <- fits.use.inh[[m]]
  ss <- extract(mod,pars = c("s_grp")) |>
    do.call(cbind, args= _) |>
    data.table()
  l <- ncol(ss)
  setnames(ss, c("V1","V2"), inhib_types, skip_absent = TRUE)
  ss <-   melt(data = ss, measure.vars=1:l)
  ss <- ss[ , .(R = quantile(value,qtl_vec),
                qntl = paste0("Q",round(100*qtl_vec,0)),
                Analyte_lab = m), by = variable]
  ss <- dcast(ss, ...~qntl, value.var = "R")
},simplify = FALSE)
smpl.use.inh <- do.call(rbind, smpl.use.inh)
smpl.use.inh[ , variable := factor(variable, inhib_types)
            ][ , grp := as.integer(variable)
            ][ , Analyte := str_remove(Analyte_lab,"\\d|(tot)")
            ][ , Lab := str_extract(Analyte_lab,"\\d|(tot)")
            ][ , facet_labs := ifelse(Lab == "tot",paste0("Σ(all)"),paste0("[M+",Lab,"]"))]

d.use.inh.agr <- dtl[ ,.(conc_mean = mean(Conc_raw_alt)), by = .(Analyte,Lab,Treatment,spl,grp) 
                    ][ , facet_labs := ifelse(Lab == "tot",paste0("Σ(all)"),paste0("[M+",Lab,"]"))]


p.use.inh <- ggplot(smpl.use.inh, aes(x = grp - 0.2, y = Q50)) +
  geom_point(size = 1) +
  geom_linerange(aes(ymin = Q25, ymax = Q75), linewidth = 0.75)+
  geom_linerange(aes(ymin = Q10, ymax = Q90), linewidth = 0.25)+
  geom_point(data = d.use.inh.agr, aes(x = grp +0.2,y = conc_mean),shape = 21, fill = "#FF000070",size =1)+
  ggh4x::facet_grid2(facet_labs~Analyte,scales = "free_y", independent = "y",)+
  labs( y = "Concentration in vial (µM)")+
  scale_x_continuous(breaks = 1:2, labels = c("Vehicle","SAM486A"), limits = c(0.5,2.5))+
  theme_bw() +
  theme(panel.grid.minor.x  = element_blank(),
        axis.line.y = element_line(linewidth  = 0.5,inherit.blank = FALSE),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
p.use.inh
ggsave("Fig08_inhibition.png",plot = p.use.inh,path =paste0(getwd(),"/graphs_tables"),width = 5, height = 6,dpi = 600 )
