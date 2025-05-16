# Combine both extraction and recovery plots
p.cmb <-
  p.ext +
  p.rec +
  plot_annotation(tag_levels = 'A')
p.cmb
ggsave("Fig04_extraction_recovery.png",plot = p.cmb,path =paste0(getwd(),"/graphs_tables"),width = 4.5, height = 4,dpi = 600 )

p.cmb2 <- 
  p.var02 + 
  plot_spacer() +
  p.var03 +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(3,0.25, 2))
p.cmb2
ggsave("Fig07_var_delta.png",plot = p.cmb2,path =paste0(getwd(),"/graphs_tables"),width = 7, height = 4,dpi = 600 )

p.qd <- 
  free(p.qnt02) +
  plot_spacer() +
  p.dev +
  plot_annotation(tag_levels = 'A') +
  plot_layout(nrow = 3, heights = c(1,0.25,1))
p.qd
ggsave("Fig05_curves_dev.png",plot = p.qd,path =paste0(getwd(),"/graphs_tables"),width = 8, height = 7,dpi = 600 )

# combine data for assay metrics table
frmt <-  function(d) signif(d, digits = 2) |> sprintf("%.12g", x=_) 
tab.mtrx <- rbind(smpl.rec,smpl.isp, fill = TRUE)
tab.mtrx[ , Assay := ifelse(grepl("R_", variable),"Recovery","Attenuation")
       ][Assay == "Recovery"     , tab_val := paste0(frmt(Q50)," (",frmt(Q25),"-",frmt(Q75),")" )
       ][Assay == "Attenuation"  , tab_val := paste0(frmt(Q50)," (",frmt(Q25),"-",frmt(Q75),")" )]
dcast(tab.mtrx, Analyte~variable, value.var = "tab_val") |>
  write.table(x = _,file = paste0(getwd(),"/graphs_tables/Table04_Metrics.csv"),sep = ",",row.names = FALSE)

# 
tab.vars <- rbind(smpl.var,smpl.mde, fill = TRUE)
tab.vars[, tab_val := paste0(frmt(Q50)," (",frmt(Q25),"-",frmt(Q75),")" )  ]
dcast(tab.vars, Analyte~metric, value.var = "tab_val") |>
  write.table(x = _,file = paste0(getwd(),"/graphs_tables/Table05_variability.csv"),sep = ",",row.names = FALSE)

