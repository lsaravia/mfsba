rm(list = ls())

# change to the folder wich contains mfSBA
#
#setwd("/cpp/mfsba")

source('Fun_MFA.r')

dq <- calcDq_mfSBA("rnd256.sed","q21.sed 2 512 20 S",F)

dq$Site <- "Random" 

hist(dq$R.Dq)

dq1<- calcDq_mfSBA("b4-991008bio.sed","q21.sed 2 512 20 S",F)

hist(dq1$R.Dq)

dq1$Site <- "Biomass"

dq <- rbind(dq,dq1)

require(ggplot2)

gp <- ggplot(dq, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp

# A publication quality graph requires lattice & Hmisc packages
#
plot_DqCI(dq)

# Plot for the spatial distribution of a spatial neutral model with 64 species 
#
plot_sed_image("t64-0100.sed","Neutral Multi species 2D distribution")

# Calculate the standardad multifractal distributions as if species represent 
# biomass or densities
#
dq1<- calcDq_mfSBA("t64-0100.sed","q21.sed 2 512 20 S",T)

dq1$Site <- "Untransformed"


dq<- calcDq_multiSBA("t64-0100.sed","q21.sed 2 512 20 S",T)
#
# Compare output with the file "testOut/s.t64-0100_SRS.sed"
#

dq$Site <- "Species Rank Surface"

dq <- rbind(dq,dq1)

gp <- ggplot(dq, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp


# Multifractal SAD
# from a neutral spatial distribution
#
dq1<- calcDq_multiSBA("t64-0100.sed","q21.sed 2 512 20 E",T)
#
# Compare output with the file "testOut/s.t64-0100_SAD.sed"
#

hist(dq1$R.Dq)

dq1$Site <- "Neutral SAD"
dq <- dq1

# From a random uniform spatial distribution of 64 species
#
dq1<- calcDq_multiSBA("t064-0256.sed","q21.sed 2 512 20 E",T)

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

# Plot of the multifractal fit for several q for the neutral spatial distribution
#
#png("Fig2.png", width=6,height=6,units="in",res=600)
plotDqFit("t.t64-0100.sed","q21.sed")
#dev.off()

# The fit of the random distribution is not good, after analysing this graph 
# we could determine that the range of boxes used should be restricted to 2-16
#
plotDqFit("t.t064-0256.sed","q21.sed")

# Estimate with different range

dq1<- calcDq_multiSBA("t064-0256.sed","q21.sed 2 16 20 E",T)
plotDqFit("t.t064-0256.sed","q21.sed")
dq1$Site <- "Rand uniform SAD"

dq <- dq[dq$Site=="Neutral SAD",]

dq <- rbind(dq,dq1)

gp <- ggplot(dq, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp

#
# For the most dominant species the distribution is more similar to a random
# form the rare species the distribution is very different from random
#

# Test Using the BCI dataset, not public, but interesting to 
# compare with paper of Borda-de-Agua AmNat 2002
# 
BCIwd<-"~/Dropbox/Projects/BCI"   
oldcd <- getwd()
setwd(BCIwd)

# Convert to sed with a resolution of 0.5 m
# 
system("./XY2sed BCI1982Sp.txt BCI1982Full.sed 1000 .5 500 .5 2")


dq1<- calcDq_multiSBA("BCI1982Full.sed","q27.sed 4 512 20 E",T)

dq1$Site <- "BCI1982 4-512"

hist(dq1$R.Dq)

gp <- ggplot(dq1, aes(x=q, y=Dq, color=Site)) +
  geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
  geom_line() +
  geom_point()
gp

plotDqFit("t.BCI1982Full.sed","q27.sed")

dq1

setwd(oldcd)

# A test for a multifractal SAD distribution with fixed Dq=0.5
# 

dq1<- calcDq_multiSBA("t08-032.sed","q21.sed 2 512 20 E",T)

dq1

plotDqFit("t.t08-032.sed","q21.sed")

# If we calculate SRS we can view the difference


dq1<- calcDq_multiSBA("t08-032.sed","q21.sed 2 512 20 S",T)

dq1

plotDqFit("t.t08-032.sed","q21.sed")

# and as the species have all the same density the standard MFA is 
# exactly the same that SRS

dq1<- calcDq_mfSBA("t08-032.sed","q21.sed 2 512 20 S",T)

dq1

# System command to delete output files
# rm a.*.sed f.*.sed s.*.sed t.*.sed 
#