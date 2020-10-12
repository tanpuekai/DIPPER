library(scales)
library(stringr)

load("G:/kathy-cheah-ddd/resources/MATRISOME.RData")
load("G:/kathy-cheah-ddd/resources/NABA.break.down.RData")
load("G:/databases/gencode/v25/gencode.v25.annotation.gene.gtf.RData")
load("G:/kathy-cheah-ddd/resources/TF.CD.receptors.human.mouse.RData")


oydata<-read.csv("old-and-young.aug21-2019.txt",sep="\t")

head(oydata)
geneSymb<-as.character(oydata[,3])
geneSymb2<-gsub(";.*$","",geneSymb)

setdiff(geneSymb2,annot[,"gene_name"])
nonMatrisome<-setdiff(geneSymb2,MATRISOME[,1])

data0<-as.matrix(oydata[,-seq(3)])

agegrp<-"old_young"
sampleIDS<-colnames(data0)
Lpheno<-strsplit(sampleIDS,"\\.")
Ages<-tolower(sapply(Lpheno,function(x)x[2]))

data<-data0[,Ages%in%unlist(strsplit(agegrp,"_"))]
sampleIDS<-colnames(data)
Lpheno<-strsplit(sampleIDS,"\\.")

Levels<-sapply(Lpheno,function(x)x[1])
Ages<-tolower(sapply(Lpheno,function(x)x[2]))
APLRC<-sapply(Lpheno,function(x)x[3])
AFNP<-sapply(Lpheno,function(x)x[4])
AFNP[is.na(AFNP)]<-APLRC[is.na(AFNP)]
AFNP2<-AFNP
AFNP2[grep("AF",AFNP)]<-"allelse"
AFNP2[AFNP=="IAF_NP"]<-"NPoriN"
AFNP2[AFNP=="NP"]<-"NPorNi"

AFNP3<-AFNP
AFNP3[grep("AF",AFNP)]<-"others"

AP<-APLRC
AP[!AP%in%c("A","P")]<-NA
LR<-APLRC
LR[!LR%in%c("L","R")]<-NA
#######################################
source("step5.helper.funcs.R")

Lgenecat<-LNABA2[c(7,5,3,2,4,6)]
Lgenecat[["nonMatrisome"]]<-setdiff(geneSymb2,unlist(Lgenecat))

f1<-factor(AFNP,levels=c("NP", "IAF_NP", "IAF", "OAF"))
ORD1<-order(paste0(Ages,".",as.integer(f1),".",APLRC,".",Levels))

MATxy<-t(apply(data[,ORD1],2,function(x){
	sapply(Lgenecat,function(y){
		length(intersect(geneSymb2[which(x>0)],y))
		})
	}))
RATIOxy<-apply(MATxy,1,function(x)x/sum(x))


#######################################
#OUTDIR<-"DEG-old-young"
OUTDIR<-paste0("DEG-for-old-young-",Suffix)
if(!dir.exists(OUTDIR))dir.create(OUTDIR)

pdf(paste0("DEG.young.vs.old.21comparisons.",Suffix,".pdf"),width=12)
	compi<-0

	p1<-VOLCANO(Ages,"young","old")
##############

	oAF<-AFNP
	oAF[oAF!="OAF"]<-"others"
	vec1<-paste0(Ages,"_",oAF)
	p2<-VOLCANO(vec1,"young_OAF","old_OAF")
	withANITA("old AF - young AF (microarray)",p2)
##############
	nonoAF<-AFNP
	nonoAF[nonoAF!="OAF"]<-"NONoAF"
	vec1<-paste0(Ages,"_",nonoAF)
	p3<-VOLCANO(vec1,"young_NONoAF","old_NONoAF")
	withANITA("old NP - young NP (microarray)",p3)

##############
	NP<-AFNP
	NP[NP!="NP"]<-"others"
	vec1<-paste0(Ages,"_",NP)
	p4<-VOLCANO(vec1,"young_NP","old_NP")
##############
	iAF<-AFNP
	iAF[iAF!="IAF"]<-"others"
	vec1<-paste0(Ages,"_",iAF)
	p6<-VOLCANO(vec1,"young_IAF","old_IAF")
##############
##############
	iNP<-AFNP
	iNP[iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP)
	p5<-VOLCANO(vec1,"young_IAF_NP","old_IAF_NP")
##############
	iNP_ap<-AFNP
	iNP_ap[iNP_ap!="IAF_NP"|is.na(AP)]<-"others"
	vec1<-paste0(Ages,"_",iNP_ap,"_AP")
	table(vec1,AP)
	p7<-VOLCANO(vec1,"young_IAF_NP_AP","old_IAF_NP_AP")
##############
	iNP_LR<-AFNP
	iNP_LR[iNP_LR!="IAF_NP"|is.na(LR)]<-"others"
	vec1<-paste0(Ages,"_",iNP_LR,"_LR")
	table(vec1,LR)
	p8<-VOLCANO(vec1,"young_IAF_NP_LR","old_IAF_NP_LR")
##############
##############
	oAF_ap<-AFNP
	oAF_ap[oAF_ap!="OAF"|is.na(AP)]<-"others"
	vec1<-paste0(Ages,"_",oAF_ap,"_AP")
	table(vec1,AP)
	p7<-VOLCANO(vec1,"young_OAF_AP","old_OAF_AP")
##############
	oAF_LR<-AFNP
	oAF_LR[oAF_LR!="OAF"|is.na(LR)]<-"others"
	vec1<-paste0(Ages,"_",oAF_LR,"_LR")
	table(vec1,LR)
	p9<-VOLCANO(vec1,"young_OAF_LR","old_OAF_LR")
##############
	oAF_ap<-AFNP
	oAF_ap[oAF_ap!="OAF"|is.na(AP)]<-"NONoAF"
	vec1<-paste0(Ages,"_",oAF_ap,"_AP")
	table(vec1,AP)
	p8<-VOLCANO(vec1,"young_NONoAF_AP","old_NONoAF_AP")
##############
	oAF_LR<-AFNP
	oAF_LR[oAF_LR!="OAF"|is.na(LR)]<-"NONoAF"
	vec1<-paste0(Ages,"_",oAF_LR,"_LR")
	table(vec1,LR)
	p10<-VOLCANO(vec1,"young_NONoAF_LR","old_NONoAF_LR")
##############
	oAF_ap<-AFNP
	oAF_ap[oAF_ap!="OAF"|is.na(AP)]<-"NONoAF"
	vec1<-paste0(Ages,"_",oAF_ap,"_AP")
	table(vec1,AP)
	p10<-VOLCANO(vec1,"young_NONoAF_AP","old_NONoAF_AP")
##############
##############
	isOAF<-AFNP
	isOAF[isOAF!="OAF"]<-"NONoAF"
	vec1<-paste0(Ages,"_",isOAF)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_NONoAF","young_OAF")
##############
	isLower<-Levels
	isLower[isLower!="L5_S1"]<-"upper"
	vec1<-paste0(Ages,"_",Levels)
	vec2<-paste0(Ages,"_",isLower)
	table(vec1,Ages)
	table(vec2,Ages)

	p10<-VOLCANO(vec1,"young_L5_S1","young_L4_5")
	p10<-VOLCANO(vec1,"young_L5_S1","young_L3_4")
	p10<-VOLCANO(vec1,"young_L3_4","young_L4_5")
	p10<-VOLCANO(vec2,"young_L5_S1","young_upper")
##############
	LRAP2<-rep(NA,length(LR))
	LRAP2[!is.na(LR)]<-"LR"
	LRAP2[!is.na(AP)]<-"AP"
	vec1<-paste0(Ages,"_",AFNP,"_",LRAP2)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_OAF_AP","young_OAF_LR")
	p10<-VOLCANO(vec1,"young_IAF_NP_AP","young_IAF_NP_LR")

##############
	ioAF<-AFNP
	ioAF[(ioAF!="OAF"&ioAF!="IAF")|is.na(LR)]<-"others"
	vec1<-paste0(Ages,"_",ioAF,"_LR")
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_LR","young_OAF_LR")
##############
	iNP<-AFNP
	iNP[iNP!="OAF"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_NP","young_OAF")
##############
	iNP<-AFNP
	iNP[iNP!="NP"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_NP","young_IAF_NP")
##############
	vec1<-paste0(Ages,"_",Levels)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_L3_4","old_L3_4")
	p10<-VOLCANO(vec1,"young_L4_5","old_L4_5")
	p10<-VOLCANO(vec1,"young_L5_S1","old_L5_S1")

dev.off()


		


