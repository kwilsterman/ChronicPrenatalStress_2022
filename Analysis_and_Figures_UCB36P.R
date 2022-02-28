library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(scales)
setwd("~/Documents/Research/Collaborative Manuscripts/Gotlieb_UCB36p/DATA/")

# RFRP-3 & TIDA Histology  -------------------------------------------------
ihc <- read.csv(file="IHC_data.csv",
                            header=TRUE, sep=",")
ihc$Day_cohort <- as.factor(ihc$Day_cohort)
ihc$grp <- as.factor(paste(ihc$Treatment, ihc$Day_cohort))

rf_n <- lm(RFRP.3_Count ~ Treatment*Day_cohort, data = ihc)
anova(rf_n)
pairs(emmeans(rf_n, ~ Treatment*Day_cohort), adjust="BH")

rf_i <- lm(RFRP.3_Mean_Intensity ~ Treatment*Day_cohort, data = ihc)
anova(rf_i)
pairs(emmeans(rf_i, ~ Treatment*Day_cohort), adjust="BH")

rf_s <- lm(RFRP.3_Mean_Size ~ Treatment*Day_cohort, data = ihc)
anova(rf_s)
pairs(emmeans(rf_s, ~ Treatment*Day_cohort), adjust="BH")

th_n <- lm(TH_counts ~ Treatment*Day_cohort, data = ihc)
anova(th_n)
pairs(emmeans(th_n, ~ Treatment*Day_cohort), adjust="BH")

th_c <- lm(TH_contacts ~ Treatment*Day_cohort, data = ihc)
anova(th_c)
pairs(emmeans(th_c, ~ Treatment*Day_cohort), adjust="BH")

ggplot() +
  geom_boxplot(data = ihc, aes(x = ordered_grp, y = RFRP.3_Mean_Size, outlier.shape = NA)) + 
  geom_jitter(data = ihc, aes(x = ordered_grp, y = RFRP.3_Mean_Size), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = ihc, aes(x = ordered_grp, y = RFRP.3_Mean_Intensity, outlier.shape = NA)) + 
  geom_jitter(data = ihc, aes(x = ordered_grp, y = RFRP.3_Mean_Intensity), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = ihc, aes(x = ordered_grp, y = RFRP.3_Count, outlier.shape = NA)) + 
  geom_jitter(data = ihc, aes(x = ordered_grp, y = RFRP.3_Count), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = ihc, aes(x = ordered_grp, y = TH_counts, outlier.shape = NA)) + 
  geom_jitter(data = ihc, aes(x = ordered_grp, y = TH_counts), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = ihc, aes(x = ordered_grp, y = TH_contacts, outlier.shape = NA)) + 
  geom_jitter(data = ihc, aes(x = ordered_grp, y = TH_contacts), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# pituitarygene expression! ------------------------------------------------
pit <- read.csv(file="pituitary_geneexp.csv",
                        header=TRUE, sep=",")
pit$day <- as.factor(pit$day)

pit_d2r <- lm(d2r ~ trt*day, data = pit)
anova(pit_d2r)
pairs(emmeans(pit_d2r, ~trt*day), adjust="BH")

pit_prl <- lm(prl ~ trt*day, data = pit)
anova(pit_prl)
pairs(emmeans(pit_prl, ~trt*day), adjust="BH")

pit_lhb <- lm(lhb ~ trt*day, data = pit)
anova(pit_lhb)
pairs(emmeans(pit_lhb, ~trt*day), adjust="BH")

ggplot() +
  geom_boxplot(data = pit, aes(x = ordered_day, y = d2r), outlier.shape = NA) + 
  geom_jitter(data = pit, aes(x = ordered_day, y = d2r), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = pit, aes(x = ordered_day, y = lhb), outlier.shape = NA) + 
  geom_jitter(data = pit, aes(x = ordered_day, y = lhb), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = pit, aes(x = ordered_day, y = prl), outlier.shape = NA) + 
  geom_jitter(data = pit, aes(x = ordered_day, y = prl), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# pituitary hormones ---------------------------------------------------
AM_data <- read.csv(file="PRL_LH_AM_data.csv",
                        header=TRUE, sep=",")
AM_data$day <- as.factor(AM_data$day)

prl_am <- lm(PRL ~ Treatment*day, data = AM_data)
anova(prl_am)
pairs(emmeans(prl_am, ~day), adjust="BH")

lh_am <- lm(LH ~ Treatment*day, data = AM_data)
anova(lh_am)
pairs(emmeans(lh_am, ~day), adjust="BH")

AM_data$plotaxis <- paste(AM_data$ordered_day, AM_data$Treatment)
ggplot() +
  geom_boxplot(data = AM_data, aes(x = plotaxis, y = PRL), outlier.shape = NA) + 
  geom_jitter(data = AM_data, aes(x = plotaxis, y = PRL), width = 0.1) +
  ylim(0,73) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot() +
  geom_boxplot(data = AM_data, aes(x = plotaxis, y = LH), outlier.shape = NA) + 
  geom_jitter(data = AM_data, aes(x = plotaxis, y = LH), width = 0.1) +
  ylim(0,3.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



PM_data <- read.csv(file="PRL_LH_PM_data.csv",
                        header=TRUE, sep=",")

PM_data$day <- as.factor(PM_data$day)

prl_pm <- lm(PRL ~ Treatment*day, data = PM_data)
anova(prl_pm)
pairs(emmeans(prl_pm, ~day), adjust="BH")

lh_pm <- lm(LH ~ Treatment*day, data = PM_data)
anova(lh_pm)
pairs(emmeans(lh_pm, ~day), adjust="BH")

PM_data$plotaxis <- paste(PM_data$ordered_day, PM_data$Treatment)
ggplot() +
  geom_boxplot(data = PM_data, aes(x = plotaxis, y = PRL), outlier.shape = NA) + 
  geom_jitter(data = PM_data, aes(x = plotaxis, y = PRL), width = 0.1) +
  ylim(0,73) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot() +
  geom_boxplot(data = PM_data, aes(x = plotaxis, y = LH), outlier.shape = NA) + 
  geom_jitter(data = PM_data, aes(x = plotaxis, y = LH), width = 0.1) +
  ylim(0,3.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())





# steroid hormone analyses --------------------------------------------------------
AM_dam_steroids <- read.csv(file="damhormones_P4_CORT_AM.csv",
                            header=TRUE, sep=",")
AM_dam_steroids$sacDay <- as.factor(AM_dam_steroids$sacDay)

## AM steroids
p4 <- lm(AM_P4_ngml ~ Trt*sacDay, data = AM_dam_steroids)
anova(p4)
pairs(emmeans(p4, ~ Trt*sacDay), adjust="BH")

cort <- lm(AM_CORT_ngml ~ Trt*sacDay, data = AM_dam_steroids)
anova(cort)
pairs(emmeans(cort, ~ Trt*sacDay), adjust="BH")

ggplot() +
  geom_boxplot(data = AM_dam_steroids, aes(x = ordered_grp, y = AM_P4_ngml), outlier.shape = NA) + 
  geom_jitter(data = AM_dam_steroids, aes(x = ordered_grp, y = AM_P4_ngml), width = 0.1) +
  ylim(0,105) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot() +
  geom_boxplot(data = AM_dam_steroids, aes(x = ordered_grp, y = AM_CORT_ngml), outlier.shape = NA) + 
  geom_jitter(data = AM_dam_steroids, aes(x = ordered_grp, y = AM_CORT_ngml), width = 0.1) +
  ylim(0,105) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_point(data = AM_dam_steroids, aes(x = AM_CORT_ngml, y = AM_P4_ngml, colour = ordered_grp)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



## PM steroids
PM_dam_steroids <- read.csv(file="damhormones_P4_CORT_PM.csv",
                            header=TRUE, sep=",")
PM_dam_steroids$PM_sample_dat <- as.factor(PM_dam_steroids$PM_sample_dat)


p4_pm <- lm(P4_ngml ~ Trt*PM_sample_dat, data = PM_dam_steroids)
anova(p4_pm)
pairs(emmeans(p4_pm, ~ Trt*PM_sample_dat), adjust="BH")

cort_pm <- lm(CORT_ngml ~ Trt*PM_sample_dat, data = PM_dam_steroids)
anova(cort_pm)
pairs(emmeans(cort_pm, ~ Trt*PM_sample_dat), adjust="BH")

ggplot() +
  geom_boxplot(data = PM_dam_steroids, aes(x = ordered_grp, y = P4_ngml), outlier.shape = NA) + 
  geom_jitter(data = PM_dam_steroids, aes(x = ordered_grp, y = P4_ngml), width = 0.1) +
  ylim(0,105) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot() +
  geom_boxplot(data = PM_dam_steroids, aes(x = ordered_grp, y = CORT_ngml), outlier.shape = NA) + 
  geom_jitter(data = PM_dam_steroids, aes(x = ordered_grp, y = CORT_ngml), width = 0.1) +
  ylim(0,105) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cor.test(PM_dam_steroids$P4_ngml[which(PM_dam_steroids$ordered_grp == "9stress")], PM_dam_steroids$mass_dayofbleed[which(PM_dam_steroids$ordered_grp == "9stress")])


# Maternal body mass data  -------------------------------------------------
maternal_mass <- read.csv(file="MaternalMass.csv",
                        header=TRUE, sep=",")
maternal_mass$ID <- as.factor(maternal_mass$ID)
maternal_mass$Treat <- as.factor(maternal_mass$Treat)

mass <- lmer(Mass ~ Day*Treat + (1|ID), data = maternal_mass)

anova(mass)
pairs(emmeans(mass, ~ Day*Treat))

maternal_mass$plot_days <- as.factor(maternal_mass$Day)
ggplot() +
  geom_boxplot(data = maternal_mass, aes(x = plot_days, y = Mass, fill = Treat), outlier.shape = NA) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Embryo Outcomes ---------------------------------------------------------

wts <- read.csv(file="SexTrtKey.csv",
                header=TRUE, sep=",")
wts$DamID <- as.factor(wts$DamID)
wts$PlacentaID <- as.factor(wts$PlacentaID)
wts$Treatment <- as.factor(wts$Treatment)
wts$Sex <- as.factor(wts$Sex)

embryo_wt <- lmer(rescale(Weight) ~ Treatment + Sex + (1|DamID), data = wts)
qqnorm(resid(embryo_wt))
qqline(resid(embryo_wt))

anova(embryo_wt)

wts$plotaxis <- as.factor(paste(wts$Sex, wts$Treatment))
ggplot() +
  geom_boxplot(data = wts, aes(x = plotaxis, y = Weight), outlier.shape = NA) + 
  geom_jitter(data = wts, aes(x = plotaxis, y = Weight), Completewidth = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Placenta Histology  -------------------------------------------------------
plhisto <-read.csv(file="Placenta_Histology.csv",
                    header=TRUE, sep=",")
key <- read.csv(file="SexTrtKey.csv",
                header=TRUE, sep=",")
sumhisto_nokey <- plhisto %>%
  group_by(DamID, PlacentaID) %>%
  summarize(LZvoid = sum(LZ_void_area_px),
            LZtotal = sum(LZ_total_area_mm2),
            JZvoid = sum(JZ_void_area_px),
            JZtotal = sum(JZ_total_area_mm2),
            Tort = sum(LZ_JZ_junction_mm))
sumhisto_nokey$DamID <- as.factor(sumhisto_nokey$DamID)
key$DamID <- as.factor(key$DamID)

sumhisto <- merge(sumhisto_nokey,key,by=c('DamID','PlacentaID'),all.x=T)
sumhisto$Pltotal <- sumhisto$LZtotal + sumhisto$JZtotal

##Models provided below. All model output was assessed using qqnorm(resid(model)) and qqline(resid(model)) for fit and anova() for output

absLZsize <- lmer(rescale(LZtotal) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)
relLZsize <- lmer(rescale(LZtotal) ~ rescale(Pltotal) + Treatment + Weight + Sex + (1|DamID), data = sumhisto)
absLZtissue <- lmer(rescale(LZvoid) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)
relLZtissue <- lmer(rescale(LZvoid) ~ rescale(LZtotal) + Treatment + Weight + Sex + (1|DamID), data = sumhisto)

absJZsize <- lmer(rescale(JZtotal) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)
relJZsize <- lmer(rescale(JZtotal) ~ rescale(Pltotal) + Treatment + Weight + Sex + (1|DamID), data = sumhisto)
absJZtissue <- lmer(rescale(JZvoid) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)
relJZtissue <- lmer(rescale(JZvoid) ~ rescale(JZtotal) + Treatment + Weight + Sex + (1|DamID), data = sumhisto)

plsize <- lmer(rescale(Pltotal) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)

absTort <- lmer(rescale(Tort) ~ Treatment + Weight + Sex + (1|DamID), data = sumhisto)
relTort <- lmer(rescale(Tort) ~ rescale(Pltotal) + Treatment + Weight + Sex + (1|DamID), data = sumhisto)

qqnorm(resid(m))
qqline(resid(m))
anova(absLZsize)

sumhisto$plotaxis <- as.factor(paste(sumhisto$Sex, sumhisto$Treatment))
ggplot() +
  geom_boxplot(data = sumhisto, aes(x = plotaxis, y = JZtotal), outlier.shape = NA) + 
  geom_jitter(data = sumhisto, aes(x = plotaxis, y = JZtotal), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot() +
  geom_boxplot(data = sumhisto, aes(x = plotaxis, y = LZtotal), outlier.shape = NA) + 
  geom_jitter(data = sumhisto, aes(x = plotaxis, y = LZtotal), width = 0.1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Placenta gene expression ------------------------------------------------
jz_geneexp <- read.csv(file="JZ_TBP_geneexp.csv",
                 header=TRUE, sep=",")
jz_geneexp$Mom <- as.factor(jz_geneexp$Mom)
jz_geneexp$Treatment <- as.factor(jz_geneexp$Treatment)

jz_phlda2 <- lmer(rescale(Phlda2) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)
jz_plii <- lmer(rescale(PLII) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)
jz_tpbpa <- lmer(rescale(Tpbpa) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)
jz_gr <- lmer(rescale(GR) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)
jz_bHSD1 <- lmer(rescale(X11bHSD1) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)
jz_bHSD2 <- lmer(rescale(X11bHSD2) ~ Treatment + Weight + Sex + (1|Mom), data = jz_geneexp)

anova(jz_phlda2)

#lz
lz_geneexp <- read.csv(file="LZ_TBP_RPLP_geneexp.csv",
                       header=TRUE, sep=",")
lz_geneexp$Mom <- as.factor(lz_geneexp$Mom)
lz_geneexp$Treatment <- as.factor(lz_geneexp$Treatment)

lz_phlda2 <- lmer(rescale(Phlda2) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)
lz_plii <- lmer(rescale(PLII) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)
lz_tpbpa <- lmer(rescale(Tpbpa) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)
lz_gr <- lmer(rescale(GR) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)
lz_bHSD1 <- lmer(rescale(X11bHSD1) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)
lz_bHSD2 <- lmer(rescale(X11bHSD2) ~ Treatment + Weight + Sex + (1|Mom), data = lz_geneexp)

anova(lz_plii)

##Interaction plots
ggplot(rawlz, aes(x = Weight, y = Tpbpa, color = Mom)) + 
  geom_point() + 
  facet_wrap(~ Treatment) +
  scale_y_log10() + 
  scale_x_log10()

##hormones
h <-read.csv(file="hormones.csv",
             header=TRUE, sep=",")

##merge JZ or LZ genes w/ histo data
histogeneexp <- merge(sumhisto_nona, rawjz ,by=c('Mom','Ind'),all.x=T)

combF <- histogeneexp %>%
  filter(Sex.x == 'F')
combM <- histogeneexp %>%
  filter(Sex.x == 'M')

##summaries
histo_bysex <- sumhisto_nona %>%
  group_by(Mom, Sex, Treatment) %>%
  summarize(avLZvoid = mean(LZvoid),
            avLZtissue = mean(LZtissue),
            avJZvoid = mean(JZvoid),
            avJZtissue = mean(JZtissue),
            avTort = mean(Tort))
write.csv(x=histo_bysex, file="histo_bysex.csv")

histo_byMom <- sumhisto_nona %>%
  group_by(Mom, Treatment) %>%
  summarize(avLZvoid = mean(LZvoid),
            avLZtissue = mean(LZtissue),
            avJZvoid = mean(JZvoid),
            avJZtissue = mean(JZtissue),
            avTort = mean(Tort))
write.csv(x=histo_byMom, file="histo_byMom.csv")

rawjz$Sex <- as.factor(rawjz$Sex)
jz_bysex <- rawjz %>%
  group_by(Mom, Sex, Treatment) %>%
  summarize(avWt = mean(Weight),
            avlog11bHSD1 = mean(log11bHSD1),
            avlog11bHSD2 = mean(log11bHSD2),
            avlogPLII = mean(logPLII),
            avlogPhlda2 = mean(logPhlda2),
            avlogTpbpa = mean(logTpbpa),
            avlogGR = mean(logGR))
write.csv(x=jz_bysex, file="jz_bysex.csv")

jz_byMom <- rawjz %>%
  group_by(Mom, Treatment) %>%
  summarize(avWt = mean(Weight),
            avlog11bHSD1 = mean(log11bHSD1),
            avlog11bHSD2 = mean(log11bHSD2),
            avlogPLII = mean(logPLII),
            avlogPhlda2 = mean(logPhlda2),
            avlogTpbpa = mean(logTpbpa),
            avlogGR = mean(logGR))
write.csv(x=jz_byMom, file="jz_byMom.csv")


lz_bysex <- rawlz %>%
  group_by(Mom, Sex, Treatment) %>%
  summarize(avlog11bHSD1 = mean(log11bHSD1),
            avlog11bHSD2 = mean(log11bHSD2),
            avlogPLII = mean(logPLII),
            avlogPhlda2 = mean(logPhlda2),
            avlogTpbpa = mean(logTpbpa),
            avlogGR = mean(logGR))
write.csv(x=lz_bysex, file="lz_bysex.csv")

lz_byMom <- rawlz %>%
  group_by(Mom, Treatment) %>%
  summarize(avlog11bHSD1 = mean(log11bHSD1),
            avlog11bHSD2 = mean(log11bHSD2),
            avlogPLII = mean(logPLII),
            avlogPhlda2 = mean(logPhlda2),
            avlogTpbpa = mean(logTpbpa),
            avlogGR = mean(logGR))
write.csv(x=lz_byMom, file="lz_byMom.csv")

wt <-read.csv(file="Maat_Weights.csv",
              header=TRUE, sep=",")
wt$ID <- as.factor(wt$ID)
wt$Treat <- as.factor(wt$Treat)

m <- lmer((1/Day.1) ~ Treat + (1|ID), data = wt)
anova(m)
qqnorm(resid(m))
qqline(resid(m))