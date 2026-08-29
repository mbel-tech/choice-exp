suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(hms)
})

p <- file.path(PROJECT_ROOT, "Figure 4 stuff script dataset image/scan sampling B3 choice exp.xlsx")
d <- read_excel(p, sheet = "data_clean")

df_plot <- d %>%
  rename(MEDIUM = `MEDIUM F.`) %>%
  mutate(
    Total_Flow = HIGH + MEDIUM + LOW,
    t_sec = as.numeric(timestamp - as.POSIXct("1899-12-31 00:00:00", tz = "UTC"), units = "secs"),
    t_sec = t_sec - min(t_sec, na.rm = TRUE)
  )

df_mean <- df_plot %>%
  mutate(t_sec = round(t_sec / 300) * 300) %>%
  group_by(t_sec) %>%
  summarise(
    mean_flow = mean(Total_Flow, na.rm = TRUE),
    sd_flow   = sd(Total_Flow,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(t_sec %% (5 * 60) == 0)

p_7c <- ggplot(df_mean, aes(x = t_sec, y = mean_flow)) +
  annotate("rect", xmin = 5 * 60,  xmax = 25 * 60,  ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  annotate("rect", xmin = 45 * 60, xmax = 65 * 60,  ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  annotate("rect", xmin = 85 * 60, xmax = 105 * 60, ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_flow - sd_flow, ymax = mean_flow + sd_flow),
    color = "#D55E00", width = 0, linewidth = 0.8
  ) +
  geom_point(color = "#D55E00", size = 2.5) +
  labs(
    # 2026-08-21 (MV): "min" only. The caption already says these are minutes
    # into the trial and that 0 is when the fish enter the arena.
    x = "min",
    y = "N° of fish in the flow",
    caption = "N.B: 0 represents the time the fish are placed in the arena"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(df_mean$t_sec, na.rm = TRUE), by = 5 * 60),
    labels = function(x) as.character(round(x / 60))
  ) +
  scale_y_continuous(breaks = seq(0, 5, by = 1)) +
  coord_cartesian(ylim = c(0, 5)) +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_text(hjust = 0, size = 30 * 0.8, margin = margin(t = 15)),
    plot.margin   = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x  = element_text(size = 33 * 0.7, face = "bold"),
    axis.title.y  = element_text(size = 33 * 0.7, face = "bold", margin = margin(r = 20)),
    axis.text.y   = element_text(size = 33 * 0.35 * 2),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 33 * 0.35 * 1.2),
    axis.line     = element_line(color = "black", linewidth = 1.2),
    panel.grid    = element_blank()
  )

out_png <- file.path(PROJECT_ROOT, "Figure 4 stuff script dataset image/Figure4_total_flow_mean_sd.png")
out_pdf <- file.path(PROJECT_ROOT, "Figure 4 stuff script dataset image/Figure4_total_flow_mean_sd.pdf")

ggsave(filename = out_png, plot = p_7c, width = 12, height = 8, dpi = 300)
cat("Saved PNG:", out_png, "\n")

ggsave(filename = out_pdf, plot = p_7c, width = 12, height = 8, device = cairo_pdf)
cat("Saved PDF:", out_pdf, "\n")
