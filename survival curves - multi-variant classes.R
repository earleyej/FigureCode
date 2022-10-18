#scratch space for survival curves on SCD stroke GWAS hits
library(survival)
library(ggplot2)
library(ggfortify)
library(gridExtra)

outFile="C:/Users/eearley/Documents/SCD - R21/stroke manuscript/figures/"
outPath=""



#### Load and set up data ####
setwd("C:/Users/eearley/Documents/SCD - R21/Stroke analysis/TOPmed/Results/cox_regression_N1333/")
pheno <- read.table("C:/Users/eearley/Documents/SCD - R21/Stroke analysis/TOPmed/Data/ssOnly_n1333_coxAge_10PCs_16Aug2021.txt",sep="\t",header=T,stringsAsFactors = F)
#survdat <- Surv(pheno$cox.age,event=pheno$Stroke_final)
#kmFit <- survfit(survdat ~ factor(pheno$Stroke_final))
#autoplot(kmFit)

#genotypes
hits <- read.table("hits.leadSNPs.genotypes.txt",sep="\t",header=F,stringsAsFactors = F)
#header
header <- read.table("n1333_samps.txt",sep="\t",stringsAsFactors = F,header=F)
colnames(hits)<-header

#convert from 1|0 to effect|non-effect
gt <- t(hits[,10:ncol(hits)])
colnames(gt) <- paste0(hits$CHROM,";",hits$POS,";",hits$REF,";",hits$ALT,";",round(or.ci$HR,digits=1),";",or.ci$CI,";",round(or.ci$MAF,digits=2))
gt.eff<-apply(gt,MARGIN=2,function(x) ifelse(x == "0|0","Non-effect/Non-effect",
                                         ifelse(x %in% c("0|1","1|0"),"Non-effect/Effect","Effect/Effect")))
#flip the state of chr6:120572043
idx=grep("chr6",colnames(gt.eff))
gt.eff[,idx]=ifelse(gt.eff[,idx] == "Non-effect/Non-effect","Effect/Effect",
                ifelse(gt.eff[,idx] == "Effect/Effect","Non-effect/Non-effect",gt.eff[,idx]))

num.vars = apply(gt.eff, 1, function(x) sum(grepl("Effect",x)))
var.tab = data.frame("NWD" = names(num.vars),
                     "num.vars" = num.vars)
x = merge(pheno,var.tab,by="NWD")
table(x$num.vars)
table(x$num.vars,x$Stroke_final)

# require >10 individuals in a class
x$num.vars = ifelse(x$num.vars >= 4, "4+",x$num.vars)
# re-factor
x$num.vars <- factor(x$num.var, levels = c("0", "1", "2", "3", "4+"))

# survival on cases + controls for all classes (but only 4 observations for 0 mutation class)
d.coxph <- survfit(Surv(cox.age,Stroke_final) ~ num.vars, data = x)
p<-autoplot(d.coxph, conf.int = F, surv.size=1)
p <- p +
  theme(legend.title= element_text(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits=c(0,1)) +
  scale_colour_discrete("")
p

png("C:/Users/eearley/Documents/SCD - R21/stroke manuscript/fig_multi_var.png",
    res=600, width=10, height=10, units = "cm")
print(p)
dev.off()


#
cox.mod <- coxph(Surv(cox.age,Stroke_final) ~ num.vars, data = x)
sum.cox.mod=summary(cox.mod)
exp(sum.cox.mod$coefficients[,2])
#
