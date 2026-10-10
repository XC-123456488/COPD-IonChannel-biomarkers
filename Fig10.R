library(ggplot2)
library(ggrepel)
library(RColorBrewer)

workDir <- getwd()
setwd(workDir)

df <- read.csv("results.csv", header = TRUE)
colnames(df)[1] <- "drug_id"
df <- df[df$rank != "Treat_vs_Control", ]
df$rank <- as.numeric(df$rank)
df$logit <- as.numeric(df$logit)
df$prob <- as.numeric(df$prob)

vol_df <- df[complete.cases(df[, c("drug_id", "rank", "logit", "prob")]), ]
top10_drugs <- vol_df[order(vol_df$rank), ][1:10, ]
vol_df$group <- ifelse(vol_df$drug_id %in% top10_drugs$drug_id, "Top10", "Others")
top10_color <- brewer.pal(10, "Set3")

p_top10 <- ggplot(vol_df, aes(x = logit, y = prob)) +
  geom_point(data = vol_df[vol_df$group == "Others", ], color = "gray50", alpha = 0.5, size = 0.8) +
  geom_point(data = vol_df[vol_df$group == "Top10", ], aes(color = drug_id), alpha = 0.8, size = 3) +
  scale_color_manual(values = top10_color, guide = "none") +
  geom_text_repel(data = vol_df[vol_df$group == "Top10", ], aes(label = drug_id), size = 3, max.overlaps = 50) +
  labs(title = "", x = "Logit Score", y = "Probability") +
  theme_classic()

pdf(file = "Figure10A.pdf", width = 6, height = 5)
print(p_top10)
dev.off()

bar_df <- df[complete.cases(df[, c("drug_id", "logit", "prob")]), ]
top30_drugs <- bar_df[order(bar_df$rank), ][1:30, ]
top30_drugs <- top30_drugs[order(-top30_drugs$logit), ]
top30_drugs$drug_id <- factor(top30_drugs$drug_id, levels = rev(top30_drugs$drug_id))
top30_drugs$order <- 1:nrow(top30_drugs)

bar <- ggplot(data = top30_drugs, aes(x = logit, y = drug_id)) +
  geom_bar(stat = "identity", width = 0.75, aes(fill = order)) +
  scale_fill_gradientn(colours = c("#e74c3c", "#f39c12", "#2ecc71",
                                   colorRampPalette(c("#3498db", "#9b59b6"))(27)), guide = "none") +
  geom_text(aes(label = sprintf("%.3f", prob), x = logit + 0.15), size = 3, hjust = 0) +
  xlim(0, max(top30_drugs$logit) + 1) +
  labs(x = "Logit Score", y = "", title = "Top30 Drugs by Logit Score") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 11), axis.text = element_text(size = 10))

pdf(file = "Figure10B.pdf", width = 5.5, height = 5.5)
print(bar)
dev.off()