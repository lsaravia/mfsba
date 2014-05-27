rm(list = ls())

# change to the folder wich contains mfSBA
#
setwd("/cpp/mfsba")

source('Fun_MFA.r')

require(ggplot2)


# A publication quality graph requires lattice & Hmisc packages
#

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


# Multifractal SAD
#
#
dq1<- calcDq_multiSBA("t64-0100.sed","q21.sed 2 512 20 E",T)

hist(dq1$R.Dq)

dq1$Site <- "Neutral SAD"
dq <- dq1

dq1<- calcDq_multiSBA("t010-0256.sed","q21.sed 2 512 20 E",T)

dq1$Site <- "Rand uniform SAD"
hist(dq1$R.Dq)

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
plotDqFit("t.t010-0256.sed","q21.sed")

