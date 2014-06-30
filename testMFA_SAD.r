# Test multifractals SRS (Species Rank Surface) 
#      and multifractal SAD (Species abundance distribution)
#

rm(list = ls())

# change to the folder wich contains mfSBA
#
#setwd("/cpp/mfsba")

source('Fun_MFA.r')

require(ggplot2)


# A publication quality graph requires lattice & Hmisc packages
#

plot_sed_image("t64-0100.sed","Neutral Multi species 2D distribution")

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

#png("Fig2.png", width=6,height=6,units="in",res=600)
plotDqFit("t.t64-0100.sed","q21.sed")
#dev.off()
plotDqFit("t.t064-0256.sed","q21.sed")


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

