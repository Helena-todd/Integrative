fsom <- fsom_rd
pctgs <- pctgs_rd

metacluster_MFIs <- t(apply(fsom$FlowSOM$data[, colsToUse], 2, function(x){
  tapply(x, fsom$metaclustering[fsom$FlowSOM$map$mapping[,1]], median)
}))

metacluster_pctgs <- t(apply(pctgs, 1, function(x){
  tapply(x, fsom$metaclustering, sum)
}))

group = samp_rd[rownames(metacluster_pctgs), "GROUP"]

p_v_meta <- rep(NA, ncol(metacluster_pctgs))
for (i in seq_len(ncol(metacluster_pctgs))) {
  data_tmp <- data.frame(Var = metacluster_pctgs[, i],Visit = group)
  fit <- aov(Var ~ Visit, data = data_tmp)
  p_v_meta[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
  #  library(ggplot2)
  #  ggplot(data_tmp) + geom_boxplot(aes(x = Visit, y = Var)) + theme_minimal()
}

p_v_meta_adj <- p.adjust(p_v_meta, "BH")
names(p_v_meta_adj) <- new_labels

new_cluster_labels <- paste0(new_labels, " (", round(p_v_meta_adj, 3), ")")
names(new_cluster_labels) <- new_labels

# means_meta <- apply(metacluster_pctgs, 2, function(x){tapply(x, group_res$groups, mean)})
# means_norm_meta <- (means - min(means))/(max(means) - min(means))


pctgs_to_plot <- data.frame(metacluster_pctgs,
                            ID = rownames(metacluster_pctgs),
                            Group = samp_rd[rownames(metacluster_pctgs), "GROUP"])
colnames(pctgs_to_plot)[1:36] <- new_labels
pctgs_to_plot <- gather(pctgs_to_plot, Metacluster, Pctg, -c(Group, ID))

pctgs_to_plot$Group <- c("non_tolerant" = "NON",
                         "primary_tolerant" = "T1",
                         "secondary_tolerant" ="T2")[pctgs_to_plot$Group]
pctgs_to_plot$Group <- factor(pctgs_to_plot$Group, levels = c("NON", "T2", "T1"))
pctgs_to_plot$Metacluster <- factor(new_cluster_labels[pctgs_to_plot$Metacluster],
                                    levels = new_cluster_labels[order(p_v_meta)])

ggplot(pctgs_to_plot) +
  geom_jitter(aes(x = Group, y = Pctg), col = "grey", height = 0, width = 0.3, size = 0.7) +
  geom_boxplot(aes(x = Group, y = Pctg), fill = "#FFFFFF00") +
  facet_wrap(~ Metacluster, scales = "free") +
  theme_minimal() +
  theme(strip.text = element_text(size=7))







### for MFIs

pdf("mfis_pvals.pdf", width = 11.69, height = 8.27)
for(marker in names(mfis_rd)){
  print(marker)

  mat2plot <- mfis_rd[[marker]]

  group = samp_rd[rownames(mat2plot), "GROUP"]

  p_v_meta <- rep(NA, ncol(mat2plot))
  for (i in seq_len(ncol(mat2plot))) {
    data_tmp <- data.frame(Var = mat2plot[, i],Visit = group)
    fit <- aov(Var ~ Visit, data = data_tmp)
    p_v_meta[i] <- summary(fit)[[1]][["Pr(>F)"]][1]
  }

  p_v_meta_adj <- p.adjust(p_v_meta, "BH")
  names(p_v_meta_adj) <- new_labels

  new_cluster_labels <- paste0(new_labels, "\n (", round(p_v_meta_adj, 3), ")")
  names(new_cluster_labels) <- new_labels

  #means_meta <- apply(metacluster_pctgs, 2, function(x){tapply(x, group_res$groups, mean)})
  #means_norm_meta <- (means - min(means))/(max(means) - min(means))


  mfis_to_plot <- data.frame(mat2plot,
                             ID = rownames(mat2plot),
                             Group = samp_rd[rownames(mat2plot), "GROUP"],
                             Batch = factor(samp_rd[rownames(mat2plot), "DATEOFCYTOFEXPERIMENT"]))
  colnames(mfis_to_plot)[1:36] <- new_labels
  mfis_to_plot <- gather(mfis_to_plot, Metacluster, Pctg, -c(Group, ID, Batch))

  mfis_to_plot$Group <- c("non_tolerant" = "NON",
                          "primary_tolerant" = "T1",
                          "secondary_tolerant" ="T2")[mfis_to_plot$Group]
  mfis_to_plot$Group <- factor(mfis_to_plot$Group, levels = c("NON", "T2", "T1"))
  mfis_to_plot$Metacluster <- factor(new_cluster_labels[mfis_to_plot$Metacluster],
                                     levels = new_cluster_labels[order(p_v_meta)])

  p <-ggplot(mfis_to_plot) +
    geom_jitter(aes(x = Group, y = Pctg, col = Batch), height = 0, width = 0.3, size = 0.7) +
    geom_boxplot(aes(x = Group, y = Pctg), fill = "#FFFFFF00") +
    facet_wrap(~ Metacluster, scales = "free") +
    theme_minimal() +
    theme(strip.text = element_text(size=7))

  print(p)
}
dev.off()
