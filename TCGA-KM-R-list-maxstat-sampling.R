#!/usr/bin/env Rscript
library(maxstat)
library(exactRankTests)
library(OIsurv)

mainpath<-'/media/rf/Work/TCGA data'
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
	genes<-c("PNPT1","DIS3L2","ZCCHC11","ZCCHC6","ZFP36","ZFP36L1","ZFP36L2","BCL2")

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
		##quantilelow<-quantile(level,c(0.25))
		##quantilehi<-quantile(level,c(0.75))
		##write(cancername,outputfile,append=T)
		##write("25%low",outputfile,append=T)
		##write(quantilelow,outputfile,append=T)
		##write("25%high",outputfile,append=T)
		##write(quantilehi,outputfile,append=T)

##filter out empty entries
		patientgenes<-patientgenes[complete.cases(patientgenes$'X_OS'),]

##using maxstat to determine cutoff, may need to random sample if sample size is too large
		survdata <- as.data.frame.matrix(patientgenes) 
		rownum <- nrow(survdata)
		if (rownum<=800) {
		cutoff <- tryCatch(maxstat.test(Surv(X_OS,X_OS_IND) ~ geneofinterest, data=survdata, smethod="LogRank", pmethod="HL"),error=function(w){message("cannot run maxstat")})
		cutnum<-as.numeric(cutoff$estimate)
		} else {
		print("random sampling")
		survdata<-survdata[sample(rownum,800), ]
		cutoff <- tryCatch(maxstat.test(Surv(X_OS,X_OS_IND) ~ geneofinterest, data=survdata, smethod="LogRank", pmethod="HL"),error=function(w){message("cannot run maxstat")})
		cutnum<-as.numeric(cutoff$estimate)}
		print(cutnum)
##error may occur if specific cancer has low death events
		if (length(cutnum)==0) {
		break}
		
##sort to low and high groups
		genequant<-function(geneofinterest){if(geneofinterest<cutnum){("low")} else if(geneofinterest>cutnum){("high")}}
		patientgenes$geneQ<-mapply(genequant,patientgenes$geneofinterest)
		head(patientgenes,5)

		quartdata<-subset(patientgenes, geneQ!="NULL")
		quartdata$geneQ<-as.character(quartdata$geneQ)
		quartdataname<-paste(genename,"-quartdata.txt",sep="")
		write.table(quartdata, quartdataname, sep="\t")

##survival curve
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
		resultmatrix=matrix(c(genename, cancername, sdiff$chisq, pval, cutnum, rownum), nrow=1, ncol=6)
##output as csv file, column 4 is pval
		write.table(resultmatrix, outputfile, sep=",", append=T, col.names = F, row.names = F)
		cat("processed: ", genename, " in ", cancername, "\n", sep="")
	}

}

##adjust p values for multiple comparisons
csvfiles <- list.files(path=mainpath, pattern="*.csv")
for (singlecsv in csvfiles) {
	fullfilename <-paste(mainpath,"/",singlecsv,sep="")
	generesult <- read.csv(file=fullfilename, header=FALSE, sep=",")
	pvalues <- generesult$V4
	fdrs <- p.adjust(pvalues,method="BH",n=length(pvalues))
	generesult$V7 <- fdrs
	oldfilename <- gsub("-output.csv","",fullfilename)
	newfilename <- paste(oldfilename,"-adj.csv",sep="")
	write.table(generesult,newfilename,append=T, col.names = F, row.names = F, sep=",")
}
