# =============================================================================
# FULL UPDATED SCRIPT (re-ordered blocks + join/filter/plot snippet included)
#
# - Imports Excel sheets:
#     "B3"                 -> data
#     "B3 data plus treat" -> metadata
# - Cleans data:
#     creates timestamp_rel from timepoint (1 = 00:00:00, +5 min each timepoint)
#     creates timestamp_obj (POSIXct) from timestamp_rel for plotting
# - Cleans metadata:
#     recodes treatment "treat" -> "exercise choice" (before joining)
# - JOINS + FILTERS + PLOTS (your snippet, adapted to use timestamp_obj we create)
# - Creates choice_exp_offeset (timepoint == 1, keep only VIDEO_ID + timestamp)
# - Saves choice_exp_offeset to Excel and prints path
# =============================================================================

# =============================================================================
# 1) LIBRARIES
# =============================================================================
pkgs <- c("readxl", "tidyverse", "writexl")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(tidyverse)
library(writexl)

# =============================================================================
# 2) IMPORT EXCEL
# =============================================================================
excel_path <- "C:/Users/martibel/OneDrive - Norwegian University of Life Sciences/Desktop/video_analysis_snapshot_choiceexp.xlsx"

data <- read_excel(excel_path, sheet = "B3")
metadata <- read_excel(excel_path, sheet = "B3 data plus treat")

cat("\nImported:\n",
    "- data:     ", nrow(data), " rows, ", ncol(data), " cols\n",
    "- metadata: ", nrow(metadata), " rows, ", ncol(metadata), " cols\n", sep = "")

# =============================================================================
# 3A) CLEAN / STANDARDIZE BEHAVIOR (data)
#     - derive relative time from timepoint
# =============================================================================
df_behavior_clean <- data %>%
  mutate(
    VIDEO_ID  = as.character(VIDEO_ID),
    timepoint = as.integer(timepoint),
    
    # Relative time in seconds: 1->0s, 2->300s, 3->600s, ...
    timestamp_rel_seconds = (timepoint - 1L) * 5L * 60L,
    
    # Human-readable HH:MM:SS
    timestamp_rel = sprintf(
      "%02d:%02d:%02d",
      timestamp_rel_seconds %/% 3600,
      (timestamp_rel_seconds %% 3600) %/% 60,
      timestamp_rel_seconds %% 60
    ),
    
    # POSIXct timestamp for plotting (anchored to today's date)
    timestamp_obj = as.POSIXct(Sys.Date(), tz = "Europe/Oslo") + timestamp_rel_seconds,
    
    # numeric coercions
    high   = suppressWarnings(as.numeric(high)),
    medium = suppressWarnings(as.numeric(medium)),
    low    = suppressWarnings(as.numeric(low)),
    calm   = suppressWarnings(as.numeric(calm)),
    tot    = suppressWarnings(as.numeric(tot)),
    
    # total flow
    Total_Flow = high + medium + low
  ) %>%
  filter(!is.na(VIDEO_ID), nzchar(VIDEO_ID))

# =============================================================================
# 3B) CLEAN / STANDARDIZE METADATA (metadata)
#     - recode treatment "treat" -> "exercise choice" BEFORE joining
# =============================================================================
df_meta <- metadata %>%
  mutate(
    video_ID  = as.character(video_ID),
    treatment = as.character(treatment),
    treatment = case_when(
      tolower(trimws(treatment)) == "treat" ~ "exercise choice",
      TRUE ~ treatment
    )
  )

# =============================================================================
# 6) JOIN + FILTER  (SNIPPET INSERTED HERE; adapted to use timestamp_obj)
# =============================================================================
df_joined <- df_behavior_clean %>%
  left_join(df_meta, by = c("VIDEO_ID" = "video_ID"))

df_final <- df_joined %>%
  filter(tolower(treatment) == "exercise choice") %>%
  filter(!is.na(timestamp_obj))

if (nrow(df_final) == 0) {
  stop("No valid rows after filtering for 'exercise choice' with valid timestamps.")
}

# =============================================================================
# 7) PLOT  (SNIPPET INCLUDED; duplicated as provided)
# =============================================================================
ggplot(df_final,
       aes(x = timestamp_obj, y = Total_Flow, group = VIDEO_ID, color = VIDEO_ID)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Total Flow by Trial (Exercise choice)",
    subtitle = "Total Flow = high + medium + low",
    x = "Time into trial",
    y = "Total Flow",
    color = "Video ID"
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "15 mins") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(data)

# 7) PLOT -------------------------------------------------------------------
ggplot(df_final,
       aes(x = timestamp_obj, y = Total_Flow, group = VIDEO_ID, color = VIDEO_ID)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Total Flow by Trial (Exercise choice)",
    subtitle = "Total Flow = high + medium + low",
    x = "Time into trial",
    y = "Total Flow",
    color = "Video ID"
  ) +
  theme_minimal() +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "15 mins") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(df_final)

# 7A) PLOT -------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(hms)

library(dplyr)
library(hms)

df_plot <- df_final %>%
  mutate(
    t_sec = as.numeric(as_hms(timestamp_rel)),  # total seconds
    Trial = factor(Trial)
  )

# 7B) PLOT -------------------------------------------------------------------


hhmmss <- function(x) {
  x <- round(x)
  sprintf("%02d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60)
}

min_lab <- function(x) paste0(round(x / 60), " min")

ggplot(df_plot,
       aes(x = t_sec, y = Total_Flow,
           group = Trial, color = Trial)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Total Flow by Trial (Exercise choice)",
    subtitle = "Total Flow = high + medium + low",
    x = "minutes into trial",
    y = "N° of fish in the flow",
    color = "Trial",
    caption = "N.B: 0 represents the time the fish are placed in the arena"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(df_plot$t_sec, na.rm = TRUE), by = 5 * 60),
    labels = min_lab
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption  = element_text(hjust = 0, size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

# 7BB) PLOT -------------------------------------------------------------------

library(dplyr)

df_top <- df_plot %>%
  group_by(t_sec) %>%
  slice_max(Total_Flow, n = 1, with_ties = FALSE) %>%
  ungroup()

ggplot(df_plot,
       aes(x = t_sec, y = Total_Flow,
           group = Trial, color = Trial)) +
  geom_line(linewidth = 1) +

  # ONE point per time, colored by the top trial
  geom_point(
    data = df_top,
    aes(x = t_sec, y = Total_Flow, color = Trial),
    size = 3,
    inherit.aes = FALSE
  ) +

  labs(
    title = "Total Flow by Trial (Exercise choice)",
    subtitle = "Points indicate the trial with the highest flow at each time",
    x = "minutes into trial",
    y = "N° of fish in the flow",
    color = "Trial",
    caption = "N.B: 0 represents the time the fish are placed in the arena"
  ) +
  theme_minimal() +
  scale_x_continuous(
    breaks = seq(0, max(df_plot$t_sec, na.rm = TRUE), by = 5 * 60),
    labels = min_lab
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption  = element_text(hjust = 0, size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )




# 7C) PLOT -------------------------------------------------------------------

library(dplyr)

df_mean <- df_plot %>%
  mutate(t_sec = round(t_sec / 300) * 300) %>%   # snap to 5-min ticks
  group_by(t_sec) %>%
  summarise(
    mean_flow = mean(Total_Flow, na.rm = TRUE),
    sd_flow   = sd(Total_Flow, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(t_sec %% (5 * 60) == 0)                 # exactly 1 datapoint per x-axis tick

p_7c <- ggplot(df_mean, aes(x = t_sec, y = mean_flow)) +
  annotate("rect", xmin = 5 * 60,  xmax = 25 * 60,  ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  annotate("rect", xmin = 45 * 60, xmax = 65 * 60,  ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  annotate("rect", xmin = 85 * 60, xmax = 105 * 60, ymin = 0, ymax = 5.5,
           fill = "#D55E00", alpha = 0.30) +
  geom_line(color = "#D55E00", linewidth = 1) +
  geom_errorbar(
    aes(ymin = mean_flow - sd_flow,
        ymax = mean_flow + sd_flow),
    color = "#D55E00",
    width = 0,
    linewidth = 0.8
  ) +
  geom_point(color = "#D55E00", size = 2.5) +
  labs(
    x = "minutes into trial",
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

print(p_7c)

ggsave(
  filename = "Figure4_total_flow_mean_sd.png",
  plot     = p_7c,
  path     = getwd(),
  width    = 12,
  height   = 8,
  dpi      = 300
)

ggsave(
  filename = "Figure4_total_flow_mean_sd.pdf",
  plot     = p_7c,
  path     = getwd(),
  width    = 12,
  height   = 8,
  device   = cairo_pdf
)


# =============================================================================
# 8) CREATE choice_exp_offeset (timepoint == 1; keep only VIDEO_ID + timestamp)
# =============================================================================
choice_exp_offeset <- df_behavior_clean %>%
  filter(timepoint == 1L) %>%
  select(VIDEO_ID, timestamp)

# =============================================================================
# 9) SAVE choice_exp_offeset AS EXCEL + PRINT PATH
# =============================================================================
out_path <- file.path(dirname(excel_path), "choice_exp_offeset.xlsx")

write_xlsx(
  x = list(choice_exp_offeset = choice_exp_offeset),
  path = out_path
)

cat("\nchoice_exp_offeset saved: ", out_path, "\n", sep = "")
