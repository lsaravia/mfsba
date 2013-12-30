rm(list = ls())

setwd("~/Dropbox/cpp/SpatialAnalysis/mfsba")

source('Fun_MFA.r')

dq <- calcDq_mfSBA("rnd256.sed","q21.sed 2 512 20 S")

dq$Site <- "Random" 

hist(dq$R.Dq)

dq1<- calcDq_mfSBA("b4-991008bio.sed","q21.sed 2 512 20 S")

hist(dq1$R.Dq)

dq1$Site <- "Biomass"

dq <- rbind(dq,dq1)

require(ggplot2)

gp <- ggplot(dq, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp


plot_DqCI(dq)

plot_sed_image("t64-0100.sed","Neutral Multi species 2D distribution")

dq1<- calcDq_mfSBA("t64-0100.sed","q21.sed 2 512 20 S",T)

dq1$Site <- "Untransformed"

dq<- calcDq_multiSBA("t64-0100.sed","q21.sed 2 512 20 S",T)

dq$Site <- "Species Rank Surface"

dq <- rbind(dq,dq1)

gp <- ggplot(dq, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp

#png("Fig1.png", width=6,height=6,units="in",res=600)
plot_DqCI(dq)
#dev.off()

#png("Fig2.png", width=6,height=6,units="in",res=600)
plotDqFit("t.t64-0100.sed","q21.sed")
#dev.off()