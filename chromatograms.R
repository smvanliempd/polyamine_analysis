#  functions
fit_gausian <- function(x,y,mu,sigma,scale,...){
  
  f = function(p){
    d = p[3]*dnorm(x,mean=p[1],sd=p[2])
    sum((d-y)^2)
  }
  
  optim(par = c(mu,sigma,scale),fn = f)
}

#get data
x <- read.csv(file=paste0(getwd(),"/data/PAs_MPs_grads_20241220.csv"),
              header = T, stringsAsFactors = F, check.names = F) |>
  data.table() |>
  melt(data = _, id.vars = 1:7,variable.name = "Analyte")

# Mobile phase optimalization with Grad1
x.mp <- x[Optimalization == "Mobile phase"
  ][between(`RT (min)`,0.45,1) & Analyte == "Putrescine" |#, value := NA
    between(`RT (min)`,0.5,1.5)  & Analyte == "Spermidine" |#, value := NA
    between(`RT (min)`,1,1.8)      & Analyte == "Spermine"  #, value := NA
  # ][`RT (min)` < 2  
  ][ , value_smt:= zoo::rollmean(value,k=2,fill = NA), by = .(`Mobile phase`,Gradient,Analyte)
  ][ , value_smt_nrm := value_smt/max(value_smt, na.rm = TRUE), by = Analyte
  ][ , Analyte := factor(Analyte, c("Putrescine","Spermidine","Spermine"))]

ggplot(x.mp, aes(x = `RT (min)`,
              y = value_smt_nrm,
              # col = Analyte
              )) +
  geom_hline(yintercept = 1, lty = 2, linewidth = 0.1)+
  geom_line(linewidth = 0.25, show.legend = FALSE) +
  scale_color_manual(values = c("#00baff","#ff4500","#ffc500"))+
  scale_x_continuous(breaks = seq(0,2,0.5),minor_breaks = seq(0,2,0.1),guide = guide_axis(minor = TRUE))+
  ggh4x::facet_grid2(Analyte ~ `Mobile phase`,scales = "free_x", independent = "x") +
  # facet_grid(Analyte ~ `Mobile phase`,scales = "free_x", space = "free_x") +
  labs(y = "Signal (max = 1 per analyte)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
ggsave("Fig01_MP_optim.png",path =paste0(getwd(),"/graphs_tables"),width = 7, height = 5,dpi = 300 )

# Gradient optimalization with MP4
x.grd <- x[ Optimalization == 'Gradient'
        ][!between(`RT (min)`,0.6,1.2) & Analyte == "Putrescine", value := NA
        ][!between(`RT (min)`,0.9,1.8) & Analyte == "Spermidine", value := NA
        ][!between(`RT (min)`,1.2,2.3) & Analyte == "Spermine"  , value := NA
        ][`RT (min)` < 2.5
        ][ , value_smt:= zoo::rollmean(value,k=2,fill = NA), by = .(`Mobile phase`,Gradient,Analyte)
        ][ , value_smt_nrm := value_smt/max(value_smt,na.rm = TRUE), by = .(Gradient,Analyte)
        ][ , Analyte := factor(Analyte, c("Putrescine","Spermidine","Spermine"))]
setkey(x.grd, "Analyte")
x.rt <- x.grd[value_smt_nrm == 1][ , .(dRT = diff(`RT (min)`)), by = Gradient][, d := paste0("d",1:.N), by = Gradient]
# dcast(x.rt, d~Gradient, value.var = "dRT")

ggplot(x.grd, aes(x = `RT (min)`,
                  y = value_smt_nrm,
                  col = Analyte)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0,2.5,0.5),minor_breaks = seq(0,2.5,0.1),guide = guide_axis(minor = TRUE))+
  ggh4x::facet_grid2(Gradient~.,scales = "free", independent = "all") +
  # facet_grid(Gradient~Analyte,scales = "free_y") +
  labs(y = "Signal (max=1 per analyte)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))

# fit gausians on peaks
x.ft <- x[ (Optimalization %in% c("Gradient")) &
          ((between(`RT (min)`,0.45,1.15) & Analyte == "Putrescine") |
           (between(`RT (min)`,0.6,1.6)  & Analyte == "Spermidine") |
           (between(`RT (min)`,1,2.1)    & Analyte == "Spermine") )
        ][ , value_smt:= zoo::rollmean(value,k=2,fill = NA)  , by = .(`Mobile phase`,Gradient,Analyte)
        ][ , value_smt_nrm := value_smt/max(value_smt,na.rm=TRUE), by = .(Gradient, `Mobile phase`, Analyte,Optimalization,`File Name`)
        ][ , Analyte := factor(Analyte, c("Putrescine","Spermidine","Spermine"))
        ][ , Optimalization := factor(Optimalization, c("Mobile phase","Gradient"))]

x.ft.pars <- x.ft[!is.na(value_smt_nrm) , {
  p <- fit_gausian(x = `RT (min)`,
                   y = value_smt_nrm, 
                   mu = .SD[value_smt_nrm == 1,`RT (min)`], 
                   sigma = 0.1, scale = 0.05,
                   control = list(maxit=10000))$par
  list(RT = p[1], pw = p[2], scl = p[3])
},by = .(`File Name`,Optimalization,  Gradient, `Mobile phase`,Analyte)]

x.ft[x.ft.pars, on = .(`File Name`,Optimalization,  Gradient, `Mobile phase`,Analyte),
     value_sim := scl * dnorm(`RT (min)`,mean = RT,sd = pw)]

ggplot(x.ft, aes(x = `RT (min)`,grp = Analyte)) +
  geom_line(aes(y = value_sim), col = "red", linewidth = .25) +
  geom_line(aes(y = value_smt_nrm), linewidth = .25) +
  scale_y_continuous(breaks = c(0,1))+
  scale_x_continuous(breaks = seq(0,2,0.5),minor_breaks = seq(0,2,0.1),guide = guide_axis(minor = TRUE))+
  facet_grid(Gradient ~ Analyte, scales = "free_x", space = "free_x" ) +
  labs(y = "Signal (max = 1 per analyte per condition)" ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"))
ggsave("Fig02_Grad_optim.png",path =paste0(getwd(),"/graphs_tables"),width = 7, height = 5,dpi = 300 )

# Make Table 3, Resolutions
setkeyv(x.ft.pars, c("Analyte","Optimalization","Mobile phase","Gradient"))
x.ft.pars[,.(R1 = (RT[2]-RT[1])/(pw[2] + pw[1]),
             R2 = (RT[3]-RT[2])/(pw[3] + pw[2])), 
          by = .(Optimalization,  Gradient)
        ][,round(.SD,1), .SDcols = c("R1","R2"),
          by = .(Optimalization,  Gradient)] |>
  write.table(x = _,file = paste0(getwd(),"/graphs_tables/Table03_Resolutions.csv"),sep = ",",row.names = FALSE)

# Plot Carry over as function of FA content in injection liquid
x.co <- x[Optimalization == "Carry over" & 
            ((between(`RT (min)`,0.45,1.15) & Analyte == "Putrescine") |
               (between(`RT (min)`,0.6,1.3)  & Analyte == "Spermidine") |
               (between(`RT (min)`,1,1.8)    & Analyte == "Spermine") )
          ][ , value_smt:= zoo::rollmean(value,k=2,fill = NA)  , by = .(Sample,Analyte)
          ][ , value_smt_nrm := value_smt/max(value_smt,na.rm=TRUE), by = Analyte
          ][ , Analyte := factor(Analyte, c("Putrescine","Spermidine","Spermine"))
          ]

ggplot(x.co, aes(x = `RT (min)`, y = value_smt, col = Sample)) +
  geom_line(linewidth = .25)+
  scale_color_manual(values = c("#00baff","#ff4500","#ffc500","#81B622","black"))+
  scale_x_continuous(breaks = seq(0,2,0.5),minor_breaks = seq(0,2,0.1),guide = guide_axis(minor = TRUE))+
  facet_grid(Analyte~.,scales = "free_y",switch = "y")+
  labs(y = "Signal (counts)")+
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.1,"cm"),
        axis.minor.ticks.length = rel(0.5),
        plot.background = element_rect(colour = "white"),
        legend.position = "inside", legend.position.inside = c(0.8,0.8)
        )
ggsave("Fig03_carry_over_FA.png",path =paste0(getwd(),"/graphs_tables"),width = 7, height = 5,dpi = 300 )
