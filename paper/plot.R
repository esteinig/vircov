library("ggplot2")
library("cowplot")
samples <- c(rep("CAMI2 (CCHFV)", 4), rep("P1 (JEV)", 4), rep("P2 (JEV)", 4), rep("P3 (LCMV)", 4))
condition <- rep(c("Confirmed", "Other"), 8)
color <- rep(c("#e31a1c", "darkgray"), 8)
analysis <- rep(c(rep("Standard", 2), rep(c("Vircov"), 2)), 4)
data <- c(6, 288, 3, 1, 2, 79, 1, 0 , 2, 111, 1, 0, 40, 95, 1, 0)
df <- data.frame(
  samples,
  condition,
  analysis,
  data,
  color
)

ggplot(df, aes(x = analysis, y = data, fill = condition)) +
  geom_bar(stat = 'identity', position = 'stack', fill=color) + facet_grid(~ samples) + theme_cowplot(12) + xlab("") + ylab("Viruses identified")

