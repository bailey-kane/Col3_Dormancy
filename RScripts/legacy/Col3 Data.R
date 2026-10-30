# load in libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggsignif)

# define functions
se <- function(x){
  sd(x)/sqrt(length(x))
}
me <- function(x, CI=0.95){
  (qt(CI,df=length(x)-1)*sd(x))/sqrt(length(x))
}

# read csv & clean
col3data <- read.csv("C:/Users/bsier/Desktop/Col1_Col3_Timepoints/Col3_Her2_DPs.csv")
col3data$ShortName = substring(col3data$SampleName, 1, 14) %>% gsub("\\.","", .)
treatment_levels <- c("MMTV-Her2 PBS",
                      "MMTV-Her2 3 dpi",
                      "MMTV-Her2 9 dpi",
                      "MMTV-Her2 28 dpi",
                      "MMTV-Her2 60 dpi",
                      "MMTV-Her2 9 mpi",
                      "IL6KO PBS",
                      "IL6KO 3 dpi",
                      "IL6KO 9 dpi",
                      "IL6KO 28 dpi")
col3data$Treatment = factor(col3data$Treatment, treatment_levels)

col3data %>%
  group_by(SampleName) %>%
  reframe(Treatment = unique(Treatment))


# summary table
gsum <- col3data %>%
  group_by(Treatment) %>%
  summarize(meanHer2count = mean(Her2PrevalencePerThousandDetections, na.rm = TRUE),
            meanDPcount = mean(DpPrevalencePerThousandDetections, na.rm = TRUE),
            meanDPpercent = mean(DoublePosPercent),
            meDPpercent = me(DoublePosPercent),
            sdDPpercent = sd(DoublePosPercent),
            meDPcount = me(DpPrevalencePerThousandDetections),
            meHer2count = me(Her2PrevalencePerThousandDetections),
            Treatment = Treatment[1]
            ) %>%
  melt(., id.vars = colnames(.[,-c(2,3)]), #preserve all cols but her2ct and dpct
       measure.vars = c("meanHer2count",
                        "meanDPcount"),
       variable.name = "Type",
       value.name = "Count")


####
#ggplots

# double bar plot with counts
ggplot(gsum, aes(Treatment, Count, fill=Type)) +
  geom_bar(stat = "identity", pos = "dodge") +
  coord_flip() +
  theme_classic() +
  ylim(0,35) +
  ylab("Count per Thousand Cells") +
  xlab("Sample ID") +
  labs(fill = "Cell type") +
  geom_text(aes(label = round(Count,2), hjust = -0.05),
            size = 3, position = position_dodge(width = 0.9)) +
  scale_fill_manual(labels = c("Her2+ Only", "Her2+:Col3a1+"),
                    values=c("#89f336", "#E0B0FF"))

# double bar plot with counts (with <0.01 for 0 values)
ggplot(gsum, aes(Treatment, Count, fill=Type)) +
  geom_bar(stat = "identity", pos = "dodge") +
  coord_flip() +
  theme_classic() +
  ylim(0,35) +
  ylab("Count per Thousand Cells") +
  xlab("Treatment") +
  labs(fill = "Cell type") +
  geom_text(aes(label = ifelse(Count >0, as.character(round(Count,2)), "<0.01")),
            hjust = -0.05, size = 3,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(labels = c("Her2+ Only", "Her2+:Col3a1+"),
                    values=c("#89f336", "#E0B0FF"))

# single overlaid bar plot with percentages
ggplot(gsum, aes(Treatment, Count, fill=Type)) +
  geom_bar(stat = "identity", pos = "identity") +
  coord_flip() +
  theme_classic() +
  ggtitle("Cell type counts") +
  ylim(0,35) +
  ylab("Count per Thousand Cells") +
  xlab("Sample ID") +
  labs(fill = "Cell type") +
  geom_text(data=gsum[gsum$Type == "meanHer2count", ],
            mapping=aes(label = paste(round(meanDPpercent,2)*100,"%",sep="")),
            hjust = -0.5, size = 3, position = position_dodge(width = 0.9)) +
  geom_text(data=gsum[gsum$Type == "meanHer2count", ],
            mapping=aes(label = paste("+/-",round(meDPpercent,2)*100,"%",sep="")),
            hjust = -1.75, size = 2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(labels = c("Her2+ Only", "Her2+:Col3a1+"),
                    values=c("#89f336", "#E0B0FF"))

#doublepos summary data
dpsum <- col3data %>%
  group_by(Treatment) %>%
  summarize(mean=mean(DoublePosPercent),
            me=me(DoublePosPercent)
            )

# bar plot with 95% CI
ggplot(dpsum, aes(Treatment, mean)) +
  geom_bar(stat = "identity", fill = "#E0B0FF") +
  theme_classic() +
  geom_errorbar(aes(x=Treatment,
                    ymin = round(mean-me,2),
                    ymax = round(mean+me,2)),
                width = 0.2) +
  geom_point(data = col3data, aes(Treatment, DoublePosPercent)) +
  coord_cartesian(ylim = c(0,1.2)) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  xlab("Sample ID") +
  ylab("Mean Proportion (95% CI)") +
  ggtitle("Double positivity as a proportion of Her2+ cells")

#box plot
ggplot(col3data, aes(Treatment, DoublePosPercent)) +
  geom_boxplot(fill = "#E0B0FF") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  xlab("Sample ID") +
  ylab("Mean Proportion") +
  ggtitle("Double positivity as a proportion of Her2+ cells")

######
#col3-ibalt proximity data

col3prox <- read.csv("D:/DORMANCY PROJECT/Exported Data/measurements_Col3.csv") %>% #reads CSV
  mutate(ShortName = substring(Image, 1, 14) %>% gsub("\\.","", .)) %>%
  mutate(Treatment = case_when( #assigns treatment to sample name (Decoding)
    ShortName == "JG-11-23_B2_6" ~ "MMTV-Her2 3 dpi",
    ShortName == "JG-22-23_C3_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-22-23_W23_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-30-23_N14_3" ~ "MMTV-Her2 9 dpi",
    ShortName == "JG-52-23_E8_3" ~ "MMTV-Her2 9 mpi",
    ShortName == "JG-59-22_F6_3" ~ "MMTV-Her2 60 dpi"
  )) %>%
mutate(Treatment = factor(Treatment, treatment_levels)) %>% #assigns factor levels to treatments so grouped by condition and time in this case
  group_by(Image) %>%
  mutate(NormCellCy5Mean = (Cell..Cy5.mean - mean(Cell..Cy5.mean, na.rm = TRUE)) / sd(Cell..Cy5.mean, na.rm = TRUE),
         NormCellTRITCMean = (Cell..TRITC.mean - mean(Cell..TRITC.mean, na.rm = TRUE)) / sd(Cell..TRITC.mean, na.rm = TRUE)) %>%
  ungroup() 


ggplot(col3prox[col3prox$Classification == "Her2+", ], aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, color = Treatment, group = Treatment)) +
    stat_smooth(data = col3prox[col3prox$Classification == "Her2-", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, color = Treatment, linetype = Classification),
              method = 'lm') +
  stat_smooth(method = "lm", aes(linetype = Classification)) +
  theme_classic() + 
  ylab("Normalized Mean Cy5 Intensity") +
  xlab("Distance to iBALT (μm)")

## col5-ibalt prox data

col5prox <- read.csv("D:/DORMANCY PROJECT/Exported Data/measurements_Col5.csv") %>% #reads CSV
  mutate(ShortName = substring(Image, 1, 14) %>% gsub("\\.","", .)) %>%
  mutate(Treatment = case_when( #assigns treatment to sample name (Decoding)
    ShortName == "JG-11-23_B2_6" ~ "MMTV-Her2 3 dpi",
    ShortName == "JG-22-23_C3_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-22-23_W23_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-30-23_N14_3" ~ "MMTV-Her2 9 dpi",
    ShortName == "JG-52-23_E8_3" ~ "MMTV-Her2 9 mpi",
    ShortName == "JG-59-22_F6_3" ~ "MMTV-Her2 60 dpi"
  )) %>%
  mutate(Treatment = factor(Treatment, treatment_levels)) %>% #assigns factor levels to treatments so grouped by condition and time in this case
  group_by(Image) %>%
  mutate(NormCellCy5Mean = (Cell..Cy5.mean - mean(Cell..Cy5.mean, na.rm = TRUE)) / sd(Cell..Cy5.mean, na.rm = TRUE),
         NormCellTRITCMean = (Cell..TRITC.mean - mean(Cell..TRITC.mean, na.rm = TRUE)) / sd(Cell..TRITC.mean, na.rm = TRUE)) %>%
  ungroup() 


ggplot(col5prox[col5prox$Classification == "Her2+", ], 
       aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean,
           color = Treatment, group = Treatment)) +
  stat_smooth(method = "lm", aes(linetype = Classification)) +
  stat_smooth(data = col5prox[col5prox$Classification == "Her2-", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, 
                  color = Treatment, linetype = Classification),
              method = 'lm') +
  theme_classic() + 
  ylab("Normalized Mean Cy5 Intensity") +
  xlab("Distance to iBALT (μm)")

# combined plot
ggplot(col5prox[col5prox$Classification == "Her2+", ], 
       aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean)) +
  stat_smooth(method = "loess", aes(linetype = Classification, color = "Col5a1")) +
  stat_smooth(data = col5prox[col5prox$Classification == "Her2-", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, 
                  color="Col5a1", linetype = Classification),
              method = 'lm') +
  stat_smooth(data = col3prox[col3prox$Classification == "Her2-", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, 
                  color="Col3a1", linetype = Classification),
              method = 'lm') +
  stat_smooth(data = col3prox[col3prox$Classification == "Her2+", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean, 
                  color = "Col3a1", linetype = Classification),
              method = 'loess') +
  theme_minimal() + 
  scale_color_manual(name = "Collagen",
                     values = c("Col3a1" = "green3",
                                "Col5a1" = "blue")) +
  facet_wrap(vars(Treatment), nrow=1) +
  ylab("Normalized Mean Cy5 Intensity (SE)") +
  xlab("Distance to iBALT (μm)")

ggplot(col3prox[sample(nrow(col3prox), 1000), ], aes(x=NormCellCy5Mean, y=NormCellTRITCMean, 
                                                     shape = Classification, color = Distance.to.annotation.iBALT.µm)) +
  geom_point() +
  coord_cartesian(ylim = c(-2,10)) +
  geom_smooth(method = "lm", inherit.aes = FALSE,
              aes(x=NormCellCy5Mean, y=NormCellTRITCMean)) +
  scale_fill_viridis_b()

# density plot of her2 cell count at each timepoint by distance
ggplot(col3prox[col3prox$Classification == "Her2+", ], aes(Distance.to.annotation.iBALT.µm, group = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.3)+
  coord_cartesian(ylim=c(0,0.006))+
  facet_wrap(~Treatment, nrow=1) +
  theme_minimal() +
  xlab("Disance to iBALT (μm)") +
  ylab("Her2+ Cell Density") +
  geom_text(col3prox[col3prox$Classification == "Her2+", ], aes(label = ))

ggplot(col5prox[col5prox$Classification == "Her2+", ], aes(Distance.to.annotation.iBALT.µm, group = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.3)+
  coord_cartesian(ylim=c(0,0.006))+
  facet_wrap(~Treatment, nrow=1) +
  theme_minimal() +
  xlab("Disance to iBALT (μm)") +
  ylab("Her2+ Cell Density")

## col3dps again
col3dat <- read.csv("D:/DORMANCY PROJECT/Exported Data/batch2_dp_counts.csv") %>%
  filter(grepl("ROI.*",Name)) %>%
  mutate(ShortName = substring(Image, 1, 14) %>% gsub("\\.","", .)) %>%
  mutate(Treatment = case_when( #assigns treatment to sample name (Decoding)
    ShortName == "JG-11-23_B2_6" ~ "MMTV-Her2 3 dpi",
    ShortName == "JG-22-23_C3_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-22-23_W23_1" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-30-23_N14_3" ~ "MMTV-Her2 9 dpi",
    ShortName == "JG-22-23_D4_7" ~ "MMTV-Her2 28 dpi",
    ShortName == "JG-52-23_E8_3" ~ "MMTV-Her2 9 mpi",
    ShortName == "JG-59-22_F6_3" ~ "MMTV-Her2 60 dpi",
    ShortName == "JG-59-22_A1_2" ~ "MMTV-Her2 PBS",
    ShortName == "JG-22_23_T20_3" ~ "IL6KO PBS",
    ShortName == "JG-12-24_B2_6" ~ "IL6KO 3 dpi",
    ShortName == "JG-12-24_I9_3" ~ "IL6KO 9 dpi",
    Image == 'JG-12-24_C3_6.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO PBS',
    Image == 'JG-12-24_E5_7.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 3 dpi',
    Image == 'JG-12-24_F6_4.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 3 dpi',
    Image == 'JG-22-23_D4_7.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 28dpi',
    Image == 'JG-22-23_U21_2.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO PBS',
    Image == 'JG-22-23_V22_1.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 28 dpi',
    Image == 'JG-23-22_A1_6.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 3 dpi',
    Image == 'JG-23-22_B2_5.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 3 dpi',
    Image == 'JG-30-23_I9_10.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 PBS',
    Image == 'JG-30-23_K11_9.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 9 dpi',
    Image == 'JG-30-23_M13_3.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 9 dpi',
    Image == 'JG-52-23_G10_9.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 9 mpi',
    Image == 'JG-52-23_H11_9.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 9 mpi',
    Image == 'JG-54-22_A1_2.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 PBS',
    Image == 'JG-54-22_B2_9.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 PBS',
    Image == 'JG-54-23_K11_7.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 9 dpi',
    Image == 'JG-59-22_H8_2.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 60 dpi',
    Image == 'JG-59-22_J10_8.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 60 dpi',
    Image == 'JG-65-22_R18_5.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 28 dpi',
    Image == 'JG-65-22_T20_6.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'IL6KO 28 dpi',
    Image == 'JG-54-23_M13_4.vsi - 20x_DAPI, TRITC, Cy5_01' ~ 'MMTV-Her2 9 dpi',
  )) %>%
  mutate(Treatment = factor(Treatment, treatment_levels),
         Her2PrevalencePerThousandDetections = (((Num.Her2...Col3a1.+Num.Her2...Col3a1..1)/Num.Detections)*1000),
         DpPrevalencePerThousandDetections = ((Num.Her2...Col3a1./Num.Detections)*1000),
         NumHer2Total = (Num.Her2...Col3a1.+Num.Her2...Col3a1..1),
         DoublePosPercent = (Num.Her2...Col3a1./(Num.Her2...Col3a1.+Num.Her2...Col3a1..1)))

col3dat[is.na(col3dat)] <- 0

tsum <- col3dat %>%
  group_by(Treatment) %>%
  summarize(meanHer2count = mean(Her2PrevalencePerThousandDetections, na.rm = TRUE),
            meanDPcount = mean(DpPrevalencePerThousandDetections, na.rm = TRUE),
            meanDPpercent = mean(DoublePosPercent),
            meDPpercent = me(DoublePosPercent),
            sdDPpercent = sd(DoublePosPercent),
            meDPcount = me(DpPrevalencePerThousandDetections),
            meHer2count = me(Her2PrevalencePerThousandDetections),
            Treatment = Treatment[1]
  ) %>%
  melt(., id.vars = colnames(.[,-c(2,3)]), #preserve all cols but her2ct and dpct
       measure.vars = c("meanHer2count",
                        "meanDPcount"),
       variable.name = "Type",
       value.name = "Count")


# double bar plot with counts (with <0.01 for 0 values)
ggplot(tsum, aes(Treatment, Count, fill=Type)) +
  geom_bar(stat = "identity", pos = "dodge") +
  coord_flip() +
  theme_classic() +
  ylim(0,45) +
  ylab("Mean Count per Thousand Cells") +
  xlab("Treatment") +
  labs(fill = "Cell type") +
  geom_text(aes(label = ifelse(Count >0, as.character(round(Count,2)), "<0.01")),
            hjust = -0.05, size = 3,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(labels = c("Her2+ Only", "Her2+:Col3a1+"),
                    values=c("#89f336", "#E0B0FF"))

#doublepos summary data
dptsum <- col3dat %>%
  group_by(Treatment) %>%
  summarize(mean=mean(DoublePosPercent),
            me=me(DoublePosPercent)
  )

# bar plot with 95% CI
ggplot(dptsum, aes(Treatment, mean)) +
  geom_bar(stat = "identity", fill = "#E0B0FF") +
  theme_classic() +
  geom_errorbar(aes(x=Treatment,
                    ymin = round(mean-me,2),
                    ymax = round(mean+me,2)),
                width = 0.2) +
  geom_point(data = col3dat, aes(Treatment, DoublePosPercent), position = position_jitter(0.05)) +
  coord_cartesian(ylim = c(0,1.2)) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  xlab("Sample ID") +
  ylab("Mean Proportion (95% CI)") +
  ggtitle("Double positivity as a proportion of Her2+ cells")

# boxplot
ggplot(col3dat, aes(Treatment, DoublePosPercent)) +
  geom_boxplot(fill = "#E0B0FF", outliers = F) +
  geom_point(shape=17, size=1) +
  stat_summary(fun=mean, geom="point", shape=15, size=3, color="purple3") +
  theme_minimal() +
  coord_cartesian(ylim = c(0,1.2)) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
  xlab("Treatment") +
  ylab("Mean Proportion (95% CI)") +
  ggtitle("Double positivity as a proportion of Her2+ cells") +
  geom_signif(comparisons = list(c("MMTV-Her2 3 dpi","MMTV-Her2 28 dpi"),
                                 c("IL6KO 9 dpi","IL6KO 28 dpi")),
              map_signif_level = TRUE,
              y_position = c(1, 1))


#####################





moremes <- read.csv("D:/DORMANCY PROJECT/Exported Data/moremes.csv") %>% #reads CSV
  group_by(Image) %>%
  mutate(NormCellCy5Mean = (Cell..Cy5.mean - mean(Cell..Cy5.mean, na.rm = TRUE)) / sd(Cell..Cy5.mean, na.rm = TRUE),
         NormCellTRITCMean = (Cell..TRITC.mean - mean(Cell..TRITC.mean, na.rm = TRUE)) / sd(Cell..TRITC.mean, na.rm = TRUE)) %>%
  ungroup() 


ggplot(moremes[moremes$Classification == "Her2+", ] , 
       aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean)) +
  stat_smooth(method = "lm", linetype = "dotted") +
  stat_smooth(data = moremes[moremes$Classification == "Her2-", ],
              aes(x=Distance.to.annotation.iBALT.µm, y=NormCellCy5Mean), 
              method = 'lm') +
  theme_classic() + 
  ylab("Normalized Mean Cy5 Intensity") +
  xlab("Distance to iBALT (μm)") +
  facet_wrap(vars(Image))

  geom_point(data = moremes[sample(nrow(moremes), 1000), ], aes(color = Classification))

ggplot(moremes[moremes$Classification == "Her2+", ], aes(Distance.to.annotation.iBALT.µm)) +
  geom_density(alpha = 0.3)+
  coord_cartesian(ylim=c(0,0.006))+
  theme_minimal() +
  xlab("Disance to iBALT (μm)") +
  ylab("Her2+ Cell Density")


# Calculate label values: (Her2+ count / total count) * 1000 per Treatment
label_data <- col3prox %>%
  group_by(Treatment) %>%
  summarise(
    total = n(),
    her2_count = sum(Classification == "Her2+")
  ) %>%
  mutate(
    value = (her2_count / total) * 1000,
    label = paste0(round(value, 2), " DCCs per 1000"),
    x = Inf,
    y = Inf
  )

# Plot
ggplot(col3prox[col3prox$Classification == "Her2+", ], 
       aes(Distance.to.annotation.iBALT.µm, group = Treatment, fill = Treatment)) +
  geom_density(alpha = 0.3) +
  coord_cartesian(ylim = c(0, 0.006)) +
  facet_wrap(~Treatment, nrow = 1) +
  theme_minimal() +
  xlab("Distance to iBALT (μm)") +
  ylab("Her2+ Cell Density") +
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3.5)
