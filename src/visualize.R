# R version of the angle distribution plot
library(ggplot2)
library(dplyr)
library(readr)

angles <- read_tsv("output/angles.tsv", show_col_types = FALSE) %>%
  mutate(angle_deg = ((angle_deg + 180) %% 360) - 180)

angles$size_cat <- factor(angles$size_cat,
                          levels = c("Tiny", "Small", "Mid", "Large", "Bulky"))

clrs <- c(
  "Tiny"  = "#4e79a7",
  "Small" = "#59a14f",
  "Mid"   = "#f28e2b",
  "Large" = "#e15759",
  "Bulky" = "#76b7b2"
)

ggplot(angles, aes(x = angle_deg, color = size_cat)) +
  geom_density(linewidth = 1.3, adjust = 1.0) +
  scale_color_manual(values = clrs, name = "Size class") +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-150, 150, 50)) +
  labs(
    title = paste0("Helix Tripeptide Angle Distribution (N=", nrow(angles), ")"),
    x = "Inter-residue vector angle (degrees)",
    y = "Density (a.u.)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "#f0f0f0"),
    plot.background  = element_rect(fill = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#999999", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    legend.position = c(0.85, 0.8)
  )

ggsave("output/distribution.png", width = 11, height = 6, dpi = 250)
