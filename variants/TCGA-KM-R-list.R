#!/usr/bin/env Rscript

mainpath<-'/media/rf/30F439E9F439B1C8/TCGA project'
outputfile<-'output.txt'
##load all folders
folderlist<-list.dirs(mainpath)
folderlist<-folderlist[-1]
print(folderlist)
for(foldername in folderlist) {
	##folderpath<-paste0(mainpath,foldername,collapse=NULL)
	print(foldername)
	setwd(foldername)
	cancername<-gsub("^.*/","",foldername)
	print(cancername)
	##write(cancername,outputfile,append=T)

##load sequencing data
	data<-read.table("genomicMatrix",sep='\t',header=T,row.names=1,check.names=F)
	head(data,5)
##filter for genes of interest
	genes<-c("PNPT1","DIS3L2","ZCCHC11","ZCCHC6")

##calculate and designate 25 and 75 quantile expression
	for(genename in genes) {
		genename
		dgenes<-data[genename,]
		tdgenes<-t(dgenes)
		head(tdgenes,5)
		colnames(tdgenes)<-"geneofinterest"

##load patient data, and merge with seq
		patient<-read.table("clinical_data",sep="\t",header=T,row.names=1)
		patientgenes<-merge(tdgenes,patient,by.x=0,by.y=0)
		head(patientgenes,5)
##setup specific gene output
		outputfile<-paste(mainpath,"/",genename,"-output.csv", sep="")
		level=patientgenes$geneofinterest
		quantilelow<-quantile(level,c(0.25))
		quantilehi<-quantile(level,c(0.75))
		##write(cancername,outputfile,append=T)
		##write("25%low",outputfile,append=T)
		##write(quantilelow,outputfile,append=T)
		##write("25%high",outputfile,append=T)
		##write(quantilehi,outputfile,append=T)
		genequart<-function(geneofinterest){if(geneofinterest<=quantilelow){("low")} else if(geneofinterest>=quantilehi){("high")}}
		patientgenes$geneQ<-mapply(genequart,patientgenes$geneofinterest)
		head(patientgenes,5)

##filter out median expression and empty entries
		quartdata<-subset(patientgenes, geneQ!="NULL")
		quartdata$geneQ<-as.character(quartdata$geneQ)
		quartdata[complete.cases(quartdata$'X_OS'),]
		head(quartdata)
		quartdataname<-paste(genename,"-quartdata.txt",sep="")
		write.table(quartdata, quartdataname, sep="\t")

##survival curve
		library(OIsurv)
##quartdata <- read.table("quartdata.txt", sep="\t", header=T)
		survgene<-Surv(as.numeric(quartdata$'X_OS'),as.numeric(quartdata$'X_OS_IND'))
		fit1 <- survfit(survgene~quartdata$geneQ)
		pdf(genename)
		plot(fit1, main="Overall Survival", xlab="Days", ylab="% Survival", col=c(1,2), lty=1)
		legend(2000, 1, c("High", "Low") , col=c(1,2), lty=1)
		dev.off()
		sdiff<-survdiff(Surv(quartdata$'X_OS',quartdata$'X_OS_IND') ~ quartdata$geneQ)
		pval<- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
		##write(c(sdiff$chisq, pval), sep=",", file=outputfile,append=T)
		resultmatrix=matrix(c(genename, cancername, sdiff$chisq, pval, quantilelow, quantilehi), nrow=1, ncol=6)
##output as csv file, column 4 is pval
		write.table(resultmatrix, outputfile, sep=",", append=T, col.names = F, row.names = F)
		cat("processed: ", genename, " in ", cancername, "\n", sep="")
	}

}
