#The heat map for the metabolite coefficients: Figure 4 from the paper Gomez et al., 2021
#Original code: Rory Wilson
#Date: 02.02.2021

rm(list=ls())

require(ggplot2)
library(dplyr)

#########################################################
#now the methylation-metabolite data

path  <- "//scidom.de/nas/AME-IGE/stat3/projects/methylation_lpsf/pfeiffer/"
npath <- paste0(path,"fig_4_wilson_yongli/")

#reading in the results

#Important:
#For your list of CpGs (N1) to include and metabolites (N2) to include,
#all results for these pairs have to be included in the results table, even if they are not significantly associated pairs.
#For example, if you wish to include N1=24 CpGs and N2=82 traits, 
#there must be 24*82 = 1968 rows of results, ie all combinations of these pairs

res    <- read.csv(file=paste0(npath,"all_results.txt"),sep="\t",stringsAsFactors = F,header=T)


#reading in the appropriate orders of metabolites, metabolite categories, CpGs and genes 
#These are the orders from top to bottom, left to right.

#Note that for the sub-orders (metabolites and CpGs): 
#the order only need to correspond to the within group/within gene orders,
#that is, for the CpG order, the actual order of the CpGs can be anything, so long as 
#that, when subsetting by gene, the order is what you want within each gene 
#the same applies to the metabolites

cpgOrd <- scan(file=paste0(npath,"order_of_cpg_sites.txt") ,what=character())
genOrd <- scan(file=paste0(npath,"order_of_genes.txt")     ,what=character())
trtOrd <- scan(file=paste0(npath,"order_of_traits.txt")    ,what=character())
grpOrd <- scan(file=paste0(npath,"order_of_trait_cats.txt"),what=character())

#making the important variables factors in the correct order
res$CpG     <- factor(res$CpG    ,rev(cpgOrd)) #this is strange, I know
res$Met     <- factor(res$Met    ,trtOrd)
res$Gene    <- factor(res$Gene   ,genOrd)
res$grp     <- factor(res$grp    ,grpOrd)

#getting better names for the significance types
res <- mutate(res,Sig = factor(case_when(
  stat_sig=="kora_only"~"Statistically\nsignificant\nin KORA only",
  stat_sig==""~"",
  stat_sig=="replicated"~"Statistically\nsignificant\nand replicated"
),levels=c("","Statistically\nsignificant\nand replicated","Statistically\nsignificant\nin KORA only")))


#Plotting:
p <- ggplot(res, aes(Met,CpG)) + 
  
  ##blaock and white theme
  theme_bw() + 
  
  #generating the tiles
  geom_tile(aes(fill = Coef),colour = "white") + 
  
  #heatmap of coefficients
  scale_fill_gradient2(low = "red",high = "blue",mid="white",midpoint=0,name="Coef.") +
  
  #no axis titles
  labs(title="",y="",x="") +
  
  #points for the statistical significance
  geom_point(mapping=aes(shape=Sig,size=Sig),color="black") + 
  
  
  #how the significance dots will appear (nothing, dot, dot), no legend
  scale_shape_manual(values=c(NA,20,20),name="",guide=F) + 
  
  #the size of the dots
  scale_size_manual(values=c(0,3,1),name="",guide=F) + 
  
  #adjusting the appearance (axes, legend)
  theme(axis.text.x  = element_text(angle=90,vjust=0.5,size=8,hjust=0),
        axis.text.y = element_text(size=8),
        legend.key.height = unit(1.1,"cm"),
        legend.key.width = unit(0.1,"cm"),
        legend.key = element_rect(colour = "white",fill="white"),axis.ticks=element_blank()) +
  
  #adjusting margins
  theme(plot.margin=unit(c(0,0,0,0),"cm")) +
  
  #facetting by gene and metabolite group
  facet_grid(Gene ~ grp,scales="free",space="free") + 
  
  #defining how the categories look and are separated 
  theme(strip.text.x = element_text(size=6, angle=90,hjust=0),
        strip.text.y = element_text(size=6,angle=0,face="italic"),
        strip.background = element_rect(colour="black",fill="gray90",linetype = "solid"),
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

################
#Printing to file
#NOTE: This is not optimized - the journal will have standards for this. This is just a demonstration.
png(paste0(npath,"fig_4_heatmap_coefficients_v02032021.png"),width = 2000, height = 800,res=150)
print(p)
dev.off()

