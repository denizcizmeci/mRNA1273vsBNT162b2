## mRNA-1273 and BNT162b2 COVID-19
## Deniz Cizmeci, Ragon Institute
## August 9, 2021


# pre-processing data  ----------------------------------------------------

ser <- read.csv("./data/systems_serology_062321.csv")

ser$group <- ser$Pregnant..Postpartum..or.Lactating.
ser$group <- factor(ser$group, levels = c("blank", "negative", "Patient not pregnant",
                                          "pregnant 13-21w", "pregnant 22-27w", "pregnant 28-33w",
                                          "pregnant 37+ weeks", "Patient lactating"))
ser2 <- subset(ser, ser$Which.vaccine.will.the.patient.be.receiving. %in% c("BNT162b2 from Pfizer/BioNTech (mRNA)",
                                                                            "mRNA-1273 from Moderna/NIH (mRNA)"))
colnames(ser2)
sero <- subset(ser2, ser2$Pregnant..Postpartum..or.Lactating. == "Patient not pregnant")
#sero$arm <- sero$Which.vaccine.will.the.patient.be.receiving.

sero <- sero[,c(1:93)]
sero$vaccine <- factor(sero$vaccine, levels = c("BNT162b2",
                                        "mRNA-1273"))

save(sero, file = "./data/sero.RData")

col.biontech <- "#81AEC6"
col.moderna <- "#FD9464"

library(tidyverse)
long <- gather(sero, "feature", "value", -c(1:7))

library(stringr)
df.feature <- data.frame(feature = colnames(sero)[8:93], 
                         assay = str_split_fixed(colnames(sero)[8:93], "_",2)[,1],
                         antigen = str_split_fixed(colnames(sero)[8:93], "_",2)[,2])

df.feature$antigen <- factor(df.feature$antigen, levels = c("S", "RBD", "NTD", "S1",   "S2", 
                                                            "Alpha", "Beta", "Gamma", 
                                                            "D614G.S", "B1.1.7", "D614G.E484K", "D614G.K417N" ))

df.feature$assay <- factor(df.feature$assay, levels = c("ADCD", "ADCP", "ADNP", 
                                                        "ADNKA.CD107a",   "ADNKA.IFNg",  "ADNKA.MIP1b",
                                                        "IgG1", "IgG2", "IgG3",    "IgM", 
                                                        "IgA1",
                                                        "FCGR2A","FCGR2B",
                                                        "FCGR3A" ,    "FCGR3B"))
df.feature$type <- "Titer"
df.feature$type[which(df.feature$assay %in% c("ADCD", "ADCP", "ADNP","ADNKA.CD107a",   "ADNKA.IFNg",  "ADNKA.MIP1b"))]  <- "Function"
df.feature$type[which(df.feature$assay %in% c("FCGR2A","FCGR2B",
                                              "FCGR3A" ,    "FCGR3B"))]  <- "FcR"
df.feature$type <- factor(df.feature$type, levels = c("Titer", "Function", "FcR"))

long <- merge(long, df.feature, by = "feature")
save(long, file = "./data/long.RData") 


# Load data ---------------------------------------------------------------

load("./data/sero.RData")
load("./data/long.RData")

col.biontech <- "#81AEC6"
col.moderna <- "#FD9464"

g.arms <- c("BNT162b2",    "mRNA-1273")
col.arms <- c(col.biontech, col.moderna)

library("ggplot2")
library("ggpubr")
library("ggsci")
library(patchwork)
library(scales)
library(rstatix)
library(tidyverse)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library("viridis")  
library(systemsseRology)
library(DMwR)

# heatmap -----------------------------------------------------------------

my.df <- sero[order(sero$vaccine),]
row.names(my.df) <- my.df$SeromYx.ID

X <- as.matrix(my.df[,c(8:93)])
#log.ind <- which(grepl("ADCD|IgG|IgA|IgM|FcR", colnames(X)))
#X[,log.ind] <- log10(X[,log.ind])
log.ind <- which(grepl("IgG|IgA|IgM|FcR", colnames(X)))
X[,log.ind] <- log10(X[,log.ind])
X <- X[, which(!grepl("D614G|B1.1.7|Alpha|Beta|Gamma", colnames(X)))]

X <- scale(X, center = TRUE, scale = TRUE)
X <- knnImputation(X)

heat.m <- X

ann_colors = list(
  vaccine = c("BNT162b2" = col.biontech, 
          "mRNA-1273" = col.moderna) )

#row.names(heat.m) <- row.names(my.df)
myBreaks <- seq(min(heat.m, na.rm = TRUE), max(heat.m, na.rm = TRUE), length.out = 100)
myColor <- inferno(99)
pdf("Fig1A_heatmap_wovariant.pdf", height = 10, width = 15)
print(pheatmap(heat.m,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               #show_colnames = FALSE,
               show_rownames = FALSE,
               #labels_row = c(),
               annotation_row = my.df[,c("vaccine"), drop = FALSE],
               #annotation_names_row = FALSE,
               annotation_legend = FALSE,
               annotation_colors = ann_colors,
               color=myColor, breaks=myBreaks,
               fontsize_row = 8,
               fontsize_col = 8,
               #angle_col = 45,
               gaps_col = c(6, 31),
               gaps_row = c(45),
               cellheight = 8,
               cellwidth = 8
               #gaps_col = c(10,seq(20,ncol(heat),20)),
))
dev.off()


# Univariates -------------------------------------------------------------

stat.wil <- long %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ vaccine)
stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
       cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
       symbols = c("***", "**", "*", " "))

temp <- subset(long, long$antigen %in% levels(long$antigen)[1:5])
titerfcr <- subset(temp, temp$type %in% c("Titer", "FcR") )
titerfcr$assay <- factor(titerfcr$assay, levels = rev(levels(titerfcr$assay)))

df_text <- subset(stat.wil, stat.wil$q < 0.05 & stat.wil$antigen %in% levels(long$antigen)[1:5])


pdf("Figure1B_univariatesWT.pdf", width = 10, height = 3)
fig1b <- ggplot(titerfcr, aes(assay, log10(value))) + 
  geom_point( shape = 21, aes(colour = vaccine),position = position_dodge(width = 0.8))+
  geom_text(data = df_text, mapping= aes(assay, y = 4.8, label = q.star),
            hjust   = 0,
            vjust   = 0.6)+
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = vaccine),
               position = position_dodge(width = 0.8), width = 0.5)+
  scale_color_manual(breaks = g.arms, values = col.arms)+scale_fill_manual(breaks = g.arms, values = col.arms)+
  ylim(1,5)+  
  ylab("log10 MFI")+ xlab("")+ coord_flip()+
  theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8),
                     strip.background=element_rect(fill="white"),
                     aspect.ratio = 1)+
  facet_grid( ~ antigen)
dev.off()

df <- subset(long, long$type == "Function" & long$antigen == "S")

pdf("functions.pdf", height = 4, width = 7)
fig1c <- ggplot(data=df, aes(x=vaccine, y=value)) +
  geom_violin(scale = "width",aes(fill = vaccine),alpha = 0.5)+ 
  geom_point(aes(fill = vaccine), shape = 21, size = 0.8, position = position_jitterdodge(0.3))+
  scale_fill_manual(breaks = g.arms,values = col.arms)+
  theme_bw()+ xlab(" ") + 
  ylab("") +
  stat_compare_means(label = "p.signif", size = 4, hide.ns = TRUE, label.x = 1.5,label.y.npc = 0.85)+
  theme(axis.text.x = element_text( color = "black",size = 6,angle = 45, hjust = 1),
        aspect.ratio = 1,
        strip.text.x = element_text(size = 8),
        strip.background=element_rect(fill="white"),
        legend.title = element_blank(), legend.position = "none") + facet_wrap(.~assay, ncol = 6, scales = "free_y")
dev.off()


pdf("Fig1BC.pdf", height = 5, width = 8)
fig1b / fig1c
dev.off()



# variant univariates -----------------------------------------------------

temp <- subset(long, long$antigen %in% levels(long$antigen)[6:8])
titerfcr <- subset(temp, temp$type %in% c("Titer", "FcR") )
titerfcr$assay <- factor(titerfcr$assay, levels = rev(levels(titerfcr$assay)))

df_text <- subset(stat.wil, stat.wil$q < 0.05 & stat.wil$antigen %in% levels(long$antigen)[6:8])


variantuni <- ggplot(titerfcr, aes(assay, log10(value))) + 
  geom_point( shape = 21, aes(colour = vaccine),position = position_dodge(width = 0.8))+
  geom_text(data = df_text, mapping= aes(assay, y = 4.8, label = q.star),
            hjust   = 0,
            vjust   = 0.6)+
  stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., group = vaccine),
               position = position_dodge(width = 0.8), width = 0.5)+
  scale_color_manual(breaks = g.arms, values = col.arms)+scale_fill_manual(breaks = g.arms, values = col.arms)+
  ylim(1,5)+  
  ylab("log10 MFI")+ xlab("")+ coord_flip()+
  theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8),
                     strip.background=element_rect(fill="white"),
                     aspect.ratio = 1)+
  facet_grid( ~ antigen)


df <- subset(long,  long$antigen %in% levels(long$antigen)[9:12])

variantfunctions <- ggplot(data=df, aes(x=vaccine, y=value)) +
  geom_violin(scale = "width",aes(fill = vaccine),alpha = 0.5)+ 
  geom_point(aes(fill = vaccine), shape = 21, size = 0.8, position = position_jitterdodge(0.3))+
  scale_fill_manual(breaks = g.arms,values = col.arms)+
  theme_bw()+ xlab(" ") + 
  ylab("") +
  stat_compare_means(label = "p.signif", size = 4, hide.ns = TRUE, label.x = 1.5,label.y.npc = 0.85)+
  theme(axis.text.x = element_text( color = "black",size = 6,angle = 45, hjust = 1),
        aspect.ratio = 1,
        strip.text.x = element_text(size = 8),
        strip.background=element_rect(fill="white"),
        legend.title = element_blank(), legend.position = "none") + facet_grid(assay~antigen, scales = "free_y")



pdf("FigUnivariate_variant2.pdf", height = 5, width = 8)
variantfunctions
dev.off()



# figure 3 variant  cor ---------------------------------------------------

sero.long <- gather(sero[,c("IgG1_S", "vaccine", "IgG1_Alpha","IgG1_Beta","IgG1_Gamma")], "variant", "value", -c("IgG1_S", "vaccine"))

col.alpha <- "#1b9e77"
col.beta <- "#d95f02"
col.gamma <- "#7570b3"

g.var <- c("IgG1_Alpha","IgG1_Beta","IgG1_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)


igg1 <- ggplot(sero.long, aes(log10(IgG1_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1, strip.background=element_rect(fill="white"), 
                                             legend.position = "none") + ylab("Variant log10(IgG1)")+
  xlab("S log10(IgG1)")+
  scale_color_manual(breaks = g.var, values = col.var) + xlim(1.5,4.5)+ ylim(1.5,4.5)


sero.long <- gather(sero[,c("IgA1_S", "vaccine", "IgA1_Alpha","IgA1_Beta","IgA1_Gamma")], "variant", "value", -c("IgA1_S", "vaccine"))
g.var <- c("IgA1_Alpha","IgA1_Beta","IgA1_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)
igA1 <- ggplot(sero.long, aes(log10(IgA1_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1, strip.background=element_rect(fill="white"), legend.position = "none") + ylab("Variant log10(IgA1)")+
  xlab("S log10(IgA1)")+
  scale_color_manual(breaks = g.var, values = col.var) + xlim(1.5,4.5)+ ylim(1.5,4.5)


sero.long <- gather(sero[,c("IgM_S", "vaccine", "IgM_Alpha","IgM_Beta","IgM_Gamma")], "variant", "value", -c("IgM_S", "vaccine"))
g.var <- c("IgM_Alpha","IgM_Beta","IgM_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)
igm <- ggplot(sero.long, aes(log10(IgM_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1, strip.background=element_rect(fill="white"), legend.position = "none") + 
  ylab("Variant log10(IgM)")+
  ylab("S log10(IgM)")+
  scale_color_manual(breaks = g.var, values = col.var) + xlim(1.5,4.5)+ ylim(1.5,4.5)


sero.long <- gather(sero[,c("FCGR2A_S", "vaccine", "FCGR2A_Alpha","FCGR2A_Beta","FCGR2A_Gamma")], "variant", "value", -c("FCGR2A_S", "vaccine"))
g.var <- c("FCGR2A_Alpha","FCGR2A_Beta","FCGR2A_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)

fc2ah <- ggplot(sero.long, aes(log10(FCGR2A_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1,strip.background=element_rect(fill="white"),  legend.position = "none") + 
  ylab("Variant log10(FCGR2A)")+
  xlab("S log10(FCGR2A)")+
  scale_color_manual(breaks = g.var, values = col.var)+ xlim(3,4.6)+ ylim(3,4.6) #+ xlim(0,40000)+ ylim(0,40000)


sero.long <- gather(sero[,c("FCGR3A_S", "vaccine", "FCGR3A_Alpha","FCGR3A_Beta","FCGR3A_Gamma")], "variant", "value", -c("FCGR3A_S", "vaccine"))
g.var <- c("FCGR3A_Alpha","FCGR3A_Beta","FCGR3A_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)
fc3av <- ggplot(sero.long, aes(log10(FCGR3A_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1,strip.background=element_rect(fill="white"),  legend.position = "none") + 
  ylab("Variant log10(FCGR3A)")+
  xlab("S log10(FCGR3A)")+
  scale_color_manual(breaks = g.var, values = col.var)+ xlim(3,4.6)+ ylim(3,4.6) #+ xlim(0,40000)+ ylim(0,40000)

sero.long <- gather(sero[,c("FCGR3B_S", "vaccine", "FCGR3B_Alpha","FCGR3B_Beta","FCGR3B_Gamma")], "variant", "value", -c("FCGR3B_S", "vaccine"))
g.var <- c("FCGR3B_Alpha","FCGR3B_Beta","FCGR3B_Gamma")
col.var <- c(col.alpha, col.beta, col.gamma)

fc3b <- ggplot(sero.long, aes(log10(FCGR3B_S), log10(value))) + 
  geom_point(aes(color = variant), size = 0.5)+ 
  stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1, strip.background=element_rect(fill="white"),legend.position = 'none') + 
  ylab("Variant log10(FCGR3B)")+
  xlab("S log10(FCGR3B)")+
  scale_color_manual(breaks = g.var, values = col.var) + xlim(3,4.6)+ ylim(3,4.6)#+ xlim(0,40000)+ ylim(0,40000)





#igg1 + igm + igA1 + fc2ah + fc3av + fc3b


pdf("f3a_v3.pdf", width = 12, height = 5)
igg1 + igm + igA1 + fc2ah + fc3av + fc3b  + plot_layout(guides = "collect") # & theme(legend.position = 'bottom')
dev.off()

pdf("f3a_v3_legend.pdf")
ggplot(sero.long, aes(log10(FCGR3B_S), log10(value))) + 
  geom_point(aes(color = variant), size = 4)+ 
  #stat_cor(method = "pearson", aes(color = variant), size = 2)+
  theme_bw() + facet_grid(.~vaccine) + theme(aspect.ratio = 1, strip.background=element_rect(fill="white")) + 
  ylab("Variant log10(FCGR3B)")+
  xlab("S log10(FCGR3B)")+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(breaks = g.var, values = col.var) + xlim(3,4.6)+ ylim(3,4.6)#+ xlim(0,40000)+ ylim(0,40000)
dev.off()


