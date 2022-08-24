## mRNA-1273 and BNT162b2 COVID-19
## Deniz Cizmeci, Ragon Institute
## June 23, 2021


# PLSDA -------------------------------------------------------------------
library(systemsseRology)
library(ggplot2)
library(DMwR)
library(patchwork)

load("./data/sero.RData")
load("./data/long.RData")

col.biontech <- "#81AEC6"
col.moderna <- "#FD9464"
g.arms <- c("BNT162b2",    "mRNA-1273")
col.arms <- c(col.biontech, col.moderna)
my_colors <- list(group = c("BNT162b2" = col.biontech,
                            "mRNA-1273" = col.moderna))

data <- sero
X <- data[,8:93]
X <- X[, which(!grepl("Alpha|Beta|Gamma|D614G.S|D614G|B1.1.7", colnames(X)))]
log.ind <- which(grepl("IgG|IgA|IgM|FcR", colnames(X)))
X[,log.ind] <- log10(X[,log.ind])
X <- scale(X, center = TRUE, scale = TRUE)
X <- knnImputation(X)
y <- factor(data$vaccine)

df_features <- data.frame(name = colnames(X),
                          label = colnames(X))
opts_plot <- list(LV_ind = c(1,2), # which LVs to plot
                  colors = my_colors,
                  y_name = "group",
                  level = 0.75) 

# Feature selection
opts_sel <- list(threshold = 0.9, n_trials = 100, return_count = TRUE)
out <- select_repeat(X, y, selector = select_lasso, options = opts_sel)

df_count <- data.frame(features = names(out$feature_count), 
                       name = names(out$feature_count), 
                       selected = out$feature_count*100/opts_sel$n_trials,
                       mark = NA)
df_count <- df_count[which(df_count$selected > 0),]
df_count <- df_count[order(-df_count$selected),]
df_count$features <- df_features$label[match(df_count$features, df_features$name)]
df_count$features <- factor(df_count$features, levels = df_count$features)
# annotation where feature is enriched
for (ind_feat in 1:nrow(df_count)) {
  tmp_mean <- rep(NA, length = nlevels(y))
  for (ind_class in 1:nlevels(y)) {
    tmp_mean[ind_class] <- mean(X[which(y == levels(y)[ind_class]),
                                  which(colnames(X) == df_count$name[ind_feat])])
  }
  df_count$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
}
df_count$mark  <- factor(df_count$mark, levels = levels(y))

plt_bar <- ggplot(data = df_count, aes(x = features, y = selected, fill = mark)) +
  scale_fill_manual(values = my_colors$group) +
  geom_bar(stat = "identity", color = "black") +
  xlab("") + geom_hline(yintercept = opts_sel$threshold*100) +
  ylab("selected (%)") +
  labs(fill = "enriched in") +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 10, angle = 80, hjust = 1, vjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
print(plt_bar)

pdf("feature_selection_repetition_lasso90.pdf", width = 10, height = 4)
print(plt_bar)
dev.off()

set.seed(123)
opts_sel <- list(threshold = 0.90, n_trials = 100, return_count = FALSE)
sel_features <- select_repeat(X, y, selector = select_lasso, options = opts_sel)


select <- function(X, y) { return(select_repeat(X, y, selector = select_lasso, options = opts_sel)) }

method = list(select = select, 
              train = train_ropls,
              predict = predict_ropls,
              score = score_accuracy)


# opts_val <- list(n_folds = 5, rf_trials = 1, pt_trials =10) # repeat with higher numbers
# vals <- validate_repeat(X, y, method, opts_val, n_trials = 100)
# save(vals, file = "vals.RData")
# # 
# plt_val <- visualize_validate(vals, options = list(y_label = "accuracy"))
# 
# pdf("val.pdf", width = 3, height = 3)
# plt_val
# dev.off()

X_sel <- X[, sel_features]
opts_model <- list(n_LV = 2)
model <- train_ropls(X_sel, y, options = opts_model)

p.scores <- visualize_ropls_scores(model, y, options = opts_plot) +
  theme(legend.title = element_blank(), legend.position = "top",
        legend.text = element_text( size = 8))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


loadsm <- as.data.frame(model@loadingMN)
vips <- as.data.frame(model@vipVn)
loads <- cbind(loadsm, vips)
loads <- loads[order(loads$`model@vipVn`),]
loads$features <- row.names(loads)
loads$features <- factor(loads$features , levels = loads$features)

loads$group <- "mRNA-1273"
loads$group[which(loads$p1 < 0)] <- "BNT162b2"

p.loads <- ggplot(data=loads, aes(x=features, y=p1))+
  #geom_bar(stat="identity", colour = "black") +coord_flip()+
  geom_bar(stat="identity", aes(fill = group), colour = "black") +coord_flip()+
  scale_fill_manual(breaks = g.arms , values = col.arms)+
  theme_classic() + ylab("LV1 loadings") + xlab("") +theme(legend.position = "none")

# p.loads <- ggplot(data=loads, aes(x=features, y=p1))+
#   geom_segment( aes(x=features, xend=features, y=0, yend=p1, color=group), size = 1) +
#   geom_point( aes(fill=group, color = group), size=6,  shape = 21) +
#   theme_light() +
#   coord_flip() +
#   geom_hline(yintercept=0)+
#   scale_fill_manual(breaks = g.arms , values = col.arms)+
#   scale_color_manual(breaks = g.arms , values = col.arms)+
#   theme_classic() + ylab("LV1 loadings") + xlab("") +theme(legend.position = "none")


pdf("plsda_lasso90.pdf", width = 7, height = 5)
p.scores  + p.loads 
#+ plt_val
dev.off()


# co-correlate network ----------------------------------------------------

library(DMwR)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(patchwork)
library(paletteer)
library(circlize)
library(RColorBrewer)
library(Hmisc)
library(tidyverse)
library(ComplexHeatmap)
library(network)
library(stringr)

col.D614G.S <- "#E64B35B2"
col.E484K.S <- "#F39B7FB2"
col.K417N.S <- "#F39B7FB2"
col.Beta <- "#7E5431FF"
col.UKD614G.S <- "#F39B7FB2"
col.UKN501Y.S <- "#F39B7FB2"
col.Alpha <- "#F39B7FB2"
col.N501Yd69.70.S <- "#F39B7FB2"

col.RBD <- "#00A087B2"
col.E484K.RBD <- "#91D1C2B2"
col.E406Q.RBD <- "#91D1C2B2"
col.F486A.RBD <- "#91D1C2B2"
col.F490K.RBD <-"#91D1C2B2"
col.HKU1.RBD <-"#91D1C2B2" 
col.Gamma <-"#91D1C2B2"
col.N487R.RBD <-"#91D1C2B2"
col.N501Y.RBD<-"#91D1C2B2" 
col.Q493R.RBD <-"#91D1C2B2"
col.SA.RBD <-"#91D1C2B2"

col.NTD <- "#3C5488B2"
col.S1 <- "#7E6148B2"
col.S2 <- "brown"

# ### Pulling some colors from different palettes
# mypal <- c(brewer.pal(8, "Pastel2") ,brewer.pal(8, "Pastel1"))
# # j.pal <- paletteer_dynamic("cartography::red.pal", 12)[1:5]
# # lj.pal <- paletteer_dynamic("cartography::red.pal", 12)[6:12]
# j.pal <- mypal[1:5]
# lj.pal <- mypal[6:12]
# hv1.pal <- paletteer_dynamic("cartography::blue.pal", 17)
# hv2.pal <- paletteer_dynamic("cartography::green.pal", 8)
# hv3.pal <- paletteer_dynamic("cartography::purple.pal", 6)
# kv.pal <- paletteer_d("nord::aurora")[1:3]
# lv.pal <- paletteer_dynamic("cartography::pink.pal", 5)
# lv2.pal <- paletteer_dynamic("cartography::brown.pal", 5)
# lv3.pal <- paletteer_dynamic("cartography::orange.pal", 8)
# lv4.pal <- paletteer_dynamic("cartography::pastel.pal", 10)
# 
# myclrs <- c(j.pal, lj.pal, hv1.pal, hv2.pal,hv3.pal, kv.pal, lv.pal,
#             lv2.pal, lv3.pal, lv4.pal, lj.pal)
# 
# 


pred <- X

res <- rcorr(as.matrix(pred), type = "spearman")

resr <- as.data.frame(res$r)
resr$selfeat <- row.names(resr)
df.r <- subset(resr, resr$selfeat %in% sel_features)

resp <- as.data.frame(res$P)
resp$selfeat <- row.names(resp)
df.p <- subset(resp, resr$selfeat %in% sel_features)


dfr <- df.r %>% gather(feature, rho, -c(selfeat))
dfp <- df.p %>% gather(feature, p, -c(selfeat))
dfp$q <- p.adjust(dfp$p, method = "BH")


identical(dfr$selfeat, dfp$selfeat)
identical(dfr$feature, dfp$feature)

df <- cbind(dfr, dfp[,c("p","q")])
df$antigen <- as.vector(str_split_fixed(df$feature, "_", 2)[,2])

df.sub <- subset(df, df$q < 0.01 & abs(df$rho) > 0.7)

link.widths <- abs(df.sub$rho)*5
#link.widths <- ifelse(abs(df.sub$rho) > 0.5, 5, 1)

link.clr <- ifelse(df.sub$rho > 0, "#e7d4e8", "#c7eae5")

#col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#003c30","#c7eae5", "white", "#e7d4e8", "#40004b"))


net <- network(df.sub, matrix.type = "edgelist", ignore.eval = FALSE, directed = FALSE)

node.cols <- network.vertex.names(net)
fc.ind <- which(grepl("FC", node.cols))
node.cols[fc.ind] <- "#8dd3c7"
titer.ind <- which(grepl("Ig", node.cols))
node.cols[titer.ind] <- "#ffffb3"
func.ind <- which(grepl("ADCD|ADCP|ADNP|ADNK", node.cols))
node.cols[func.ind] <- "#bebada"

# node.cols <- as.vector(str_split_fixed(node.cols, "_", 2)[,2])
# node.cols <- gsub("D614G.S", col.D614G.S, node.cols)
# node.cols <- gsub("Alpha", col.Alpha, node.cols)
# node.cols <- gsub("Beta", col.Beta, node.cols)
# node.cols <- gsub("Gamma", col.Gamma, node.cols)
# node.cols <- gsub("RBD", col.RBD, node.cols)
# node.cols <- gsub("S1", col.S1, node.cols)
# node.cols <- gsub("S2", col.S2, node.cols)
# node.cols <- gsub("NTD", col.NTD, node.cols)

#node.cols[c(6, 17,28,39)] <- "#E6F5C9"

pdf("cocornet.pdf", width = 4, height = 4)
print(plot.network(net, label = network.vertex.names(net), 
                   #vertex.col = "slategray",
                   vertex.col = as.color(node.cols),
                   #mode = "circle",
                   edge.col = link.clr,
                   edge.lwd = link.widths,
                   pad = 2,
                   label.cex = 0.4))
dev.off()

