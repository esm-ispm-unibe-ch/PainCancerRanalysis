
library(devtools)
install_github("esm-ispm-unibe-ch/NMAJags")
library(NMAJags)
library(R2jags)
library(netmeta)
library(meta)
library(metafor)
library(readxl)
####




##################################################################################
#
#   DATA LOADING AND CLEANING 
#
##################################################################################


#Load the data
#make sure missing data is in empty cells!
CancerPainDATA <-read_excel("~/_mydrive/WHO/Cancer Pain/KQ 1.1-1.3 NMA pain data 1409.xlsx",na="empty")
CancerPainDATA1.2=subset(CancerPainDATA,`KQ 1.1 or 1.2`=="1.2 Maintenance")
CancerPainDATA1.2=subset(CancerPainDATA1.2,Arm!="zz_Total")
#Load the continuous data from Orestis
CancerPainDATA1.2cont <- read_excel("~/_mydrive/WHO/Cancer Pain/CONTINUOUS.imputed.data.xlsx", 
                                    col_types = c("text", "numeric", "numeric", 
                                                  "numeric", "numeric", "text", "numeric", 
                                                  "text", "text", "text", "numeric", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "numeric", "text", "numeric", 
                                                  "numeric", "text", "text", "text", 
                                                  "text", "text", "text", "text", "numeric", 
                                                  "numeric", "text", "numeric", "text", 
                                                  "text", "text", "numeric", "text", 
                                                  "numeric", "text", "text", "text", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "text", "numeric", "text", 
                                                  "numeric", "numeric", "text", "numeric", 
                                                  "text", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "text", "numeric", 
                                                  "text", "numeric", "text", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "text", "text", 
                                                  "text", "text", "text", "numeric", 
                                                  "numeric", "numeric"))

#Exclude studies with Butorphanol
CancerPainDATA1.2=subset(CancerPainDATA1.2,treatment.1!="Butorphanol")
CancerPainDATA1.2cont=subset(CancerPainDATA1.2cont,treatment.1!="Butorphanol")
#Load the treatment re-coding excel
Drug_definition <- read_excel("~/_mydrive/WHO/Cancer Pain/Drug definition.xlsx")

#recode the treatments into a new column "Treatment"
CancerPainDATA1.2$Treatment=NA
CancerPainDATA1.2cont$Treatment=NA
for (i in 1:dim(Drug_definition)[1])
  {
  CancerPainDATA1.2$Treatment[CancerPainDATA1.2$treatment.1==Drug_definition$Detailed[i]]=Drug_definition$Grouped[i]
  CancerPainDATA1.2cont$Treatment[CancerPainDATA1.2cont$treatment.1==Drug_definition$Detailed[i]]=Drug_definition$Grouped[i]
  }

#keep only useful variables

CancerPainDATA1.2dich=subset(CancerPainDATA1.2,select=c(`Study ID`,Year,`Age: Unique`, Arm,`Outcome Limit`,`OVERALL QUALITY`, Scale, Timepoint, `Timepoint Units`,`N Analyzed`,responders,Treatment))
CancerPainDATA1.2cont=subset(CancerPainDATA1.2cont,select=c(Study.ID,Year,Age..Unique, Arm,OVERALL.QUALITY, Scale, Timepoint, Timepoint.Units,N.Analyzed,responders,Mean..norm, Mean.Chg..norm,SD.change.norm,SD.norm,Treatment))


#create the continious outcome y, sd, ncont
# Get a mixture of change scores and final values
#-----------------------
  #give priority to the change data
  CancerPainDATA1.2cont$y=CancerPainDATA1.2cont$Mean.Chg..norm 
  CancerPainDATA1.2cont$ncont=CancerPainDATA1.2cont$N.Analyzed
  CancerPainDATA1.2cont$sd=CancerPainDATA1.2cont$SD.change.norm
  #if change data not reported use final values
  use.mean=is.na(CancerPainDATA1.2cont$Mean.Chg..norm)|is.na(CancerPainDATA1.2cont$sd)
  CancerPainDATA1.2cont$y[use.mean]=CancerPainDATA1.2cont$Mean..norm[use.mean]
  CancerPainDATA1.2cont$sd[use.mean]=CancerPainDATA1.2cont$SD.norm[use.mean]


#define 2 datasets, for dichotomous first************
#*** THIS IS THE OUTPUT OF THIS FIRST PART 
CancerPainDATA1.2dich=subset(CancerPainDATA1.2dich,!is.na(responders))
CancerPainDATA1.2cont=subset(CancerPainDATA1.2cont,(!is.na(sd))&(!is.na(y)))


##################################################################################
#
#   Analysis of dichotomous outcome 
#
##################################################################################

#USING THE CancerPainDATA1.2dich
#number of studies
length(table(CancerPainDATA1.2dich$`Study ID`))
#number of arms
sum(table(CancerPainDATA1.2dich$`Study ID`))
table(CancerPainDATA1.2dich$Treatment)


#######################################
#### Setup a connected network
#######################################

CancerPainDATA1.2dich=pool.arms(data=CancerPainDATA1.2dich, studlab=`Study ID`,treat=Treatment, r=responders,n=`N Analyzed`)
DATApairsd=pairwise(treat=treat,event=r,n=n, data=CancerPainDATA1.2dich, studlab=id, sm= "OR")
#get connected network
netconnection(treat1,treat2,studlab, data=DATApairsd, subset=!is.na(TE))


#######################################
#### Finally, NMA
#######################################

#run NMA and create an object called EFF for efficacy
EFFd<-netmeta(TE,seTE,treat1,treat2,studlab,data=DATApairsd,  sm="OR",r="Placebo",comb.fixed =F, comb.random = T, tol.multiarm=T,subset=!is.na(TE))
#netgraph(EFFd, plastic=F, thickness="number.of.studies", multiarm = F, points=T,cex.points=4, number.of.studies =T, col=1)
forest(EFFd, ref="placebo", sortvar = Pscore,xlab="OR")

netsplit(EFFd) 
decomp.design(EFFd)


##################################################################################
#
#   Analysis of continuous outcome 
#
##################################################################################

#USING THE CancerPainDATA1.2dich
#number of studies
length(table(CancerPainDATA1.2cont$Study.ID))
#number of arms
sum(table(CancerPainDATA1.2cont$Study.ID))
table(CancerPainDATA1.2cont$Treatment)


#######################################
#### Setup a connected network
#######################################

CancerPainDATA1.2cont=pool.arms(data=CancerPainDATA1.2cont, studlab=Study.ID,treat=Treatment, y=y,n=N.Analyzed,sd=sd, type = "cont")
DATApairsc=pairwise(treat=treat,mean=y,n=n, sd=sd, data=CancerPainDATA1.2cont, studlab=id, sm= "SMD")
#get connected network
netconnection(treat1,treat2,studlab, data=DATApairsc, subset=!is.na(TE))


#######################################
#### Finally, NMA
#######################################

#run NMA and create an object called EFF for efficacy

EFFc<-netmeta(TE,seTE,treat1,treat2,studlab,data=DATApairsc,  sm="SMD",r="Placebo",comb.fixed =F, comb.random = T, tol.multiarm=T,subset=!is.na(TE), tau.preset=sqrt(0.1))
netgraph(EFFc, plastic=F, thickness="number.of.studies", multiarm = F, points=T,cex.points=4, number.of.studies =T, col=1)
forest(EFFc, ref="Placebo",sort=-Pscore, xlab="SMD")

netsplit(EFFc) 
decomp.design(EFFc)


