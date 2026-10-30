# load in libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)

datapath <- '/Volumes/Untitled/DORMANCY PROJECT/Exported Data/measurements_Col3.csv'

# import slides decoded and exported data
decoded <- read.csv("../../Slides_Decoded_20261030.csv")
dat <- read.csv(datapath)

# use stringr to clean up slide and sample names, then add metadata from decoded
dat <- dat |>
  mutate(
    SlideName = str_remove(Image, "\\.vsi.*$"),
    SampleNumber = as.integer(str_extract(Image, "\\d+$"))
  ) |>
  left_join(
    decoded |>
      mutate(SampleNumber = as.integer(SampleNumber)) |>
      select(SlideName, SampleNumber, Model:DPI),
    by = c("SlideName", "SampleNumber")
  ) |> 
  mutate(Arm = paste(TrxGroup, DPI))

# normalize cy5 values by image
dat <- dat |> 
  group_by(Image) |> 
  mutate(NormCy5Mean = (Cell..Cy5.mean - 
                            mean(Cell..Cy5.mean, na.rm = TRUE)) 
                            / sd(Cell..Cy5.mean, na.rm = TRUE)
         ) |>
  ungroup()

# plot for multi arm
arm_order <- dat %>%
  mutate(DPI = as.integer(DPI)) %>%
  distinct(Arm, TrxGroup, DPI) %>%
  arrange(TrxGroup, DPI) %>%
  pull(Arm) %>%
  unique()

dat$Arm <- factor(dat$Arm, levels = arm_order)

arm_labels <- dat %>%
  group_by(Arm) %>%
  summarise(n = n_distinct(Image)) %>%
  mutate(label = paste0(Arm, " (n=", n, ")")) %>%
  { setNames(.$label, .$Arm) }

datrm <- dat |> filter(SlideName != "JG-22-23_C3_1")

ggplot(datrm, aes(x = Distance.to.annotation.iBALT.µm,
                y = NormCy5Mean,
                color = Classification)) +
  geom_point(data = ~ filter(.x, Classification == "Her2+"),
             color = "green1", alpha = 0.3, size = 0.25) +
  stat_smooth(data = ~ filter(.x, Classification == "Her2-"), method = "lm") +
  stat_smooth(data = ~ filter(.x, Classification == "Her2+"), method = loess) +
  facet_wrap(~ Arm, nrow = 1, labeller = labeller(Arm = arm_labels)) +
  scale_color_manual(values = c("Her2+" = "limegreen", "Her2-" = "grey60")) +
  scale_x_continuous(breaks = scales::breaks_width(100)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
        strip.background = element_rect(fill = "grey90", color = NA)) +
  ylab("Per-cell mean Cy5 intensity\n(z-score normalized)") +
  xlab("Distance to iBALT (μm)") +
  labs(title = "Col3 expression as a function of iBALT distance")
  
