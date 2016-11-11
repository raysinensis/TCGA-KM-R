#!/usr/bin/env Rscript

mainpath<-'/media/rf/Work/TCGA data/'
outputfile<-'/media/rf/Work/TCGA data/output.txt'
##load all folders
folderlist<-list.dirs(mainpath)
folderlist<-folderlist[-1]
print(folderlist)
for(foldername in folderlist) {
	##folderpath<-paste0(mainpath,foldername,collapse=NULL)
	##print(foldername)
	setwd(foldername)
	cancername<-gsub("^.*//","",foldername)
	write(cancername,outputfile,append=T)

##load sequencing data
	data<-read.table("genomicMatrix",sep='\t',header=T,row.names=1,check.names=F)
	head(data,5)
##filter for genes of interest
	genes<-c("PNPT1","DIS3L2","ZCCHC11","ZCCHC6","ZFP36","RBL2")
	dgenes<-data[genes,]
	tdgenes<-t(dgenes)
	head(tdgenes,5)

##load patient data, and merge with seq
	patient<-read.table("clinical_data",sep="\t",header=T,row.names=1)
	patientgenes<-merge(tdgenes,patient,by.x=0,by.y=0)
	head(patientgenes,5)

##calculate and designate 25 and 75 quantile expression
	level=patientgenes$PNPT1
	quantilelow<-quantile(level,c(0.25))
	quantilehi<-quantile(level,c(0.75))
	write(quantilelow,outputfile,append=T)
	write(quantilehi,outputfile,append=T)
	PNPT1quart<-function(PNPT1){if(PNPT1<=quantilelow){("low")} else if(PNPT1>=quantilehi){("high")}}
	patientgenes$PNPT1q<-mapply(PNPT1quart,patientgenes$PNPT1)
	head(patientgenes,5)

##filter out median expression and empty entries
	quartdata<-subset(patientgenes, PNPT1q!="NULL")
	quartdata$PNPT1q<-as.character(quartdata$PNPT1q)
	quartdata[complete.cases(quartdata$'X_OS'),]
	head(quartdata)
	write.table(quartdata, "quartdata.txt", sep="\t")

##survival curve
	library(OIsurv)
##quartdata <- read.table("quartdata.txt", sep="\t", header=T)
	survgene<-Surv(as.numeric(quartdata$'X_OS'),as.numeric(quartdata$'X_OS_IND'))
	fit1 <- survfit(survgene~quartdata$PNPT1q)
	pdf(cancername)
	plot(fit1, main="Overall Survival", xlab="Days", ylab="% Survival", col=c(1,2), lty=1)
	legend(2000, 1, c("High", "Low") , col=c(1,2), lty=1)
	dev.off()
	sdiff<-survdiff(Surv(quartdata$'X_OS',quartdata$'X_OS_IND') ~ quartdata$PNPT1q)
	pval<- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
	write(sdiff$chisq,outputfile,append=T)
	write(pval,outputfile,append=T)
	cat("processed: ", cancername, "\n")

}
