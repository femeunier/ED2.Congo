rm(list = ls())

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ED2.Congo)

################################################################################
# census data

data.file <- "./data/Yangambi data_by_individual main_recensus.xlsx"
data <- read_xlsx(data.file,sheet = "Data")

################################################################################
# plot data

census.file <- "./data/2_data_by_censusinterval.csv"
data.census.interval <- read.csv(census.file)

################################################################################
# Plot/census data

data.plot <- data %>%
  dplyr::select(PlotID,`Plot Code`) %>%
  distinct()
data.census <- data %>%
  dplyr::select(`Census Date`,`Census No`) %>%
  distinct() %>%
  rename(census.no = `Census No`,
         census.date = `Census Date`) %>%
  group_by(census.no) %>%
  summarise(census.date = mean(census.date),
            .groups = "keep")

################################################################################
# Individual trees

data.selected <- data %>%
  dplyr::select(TreeID,PlotID,`Census No`,`Sub Plot T1`,Species,D4) %>%
  rename(patch = `Sub Plot T1`,
         census = `Census No`,
         species = Species) %>%
  mutate(patch = as.numeric(patch)) %>%
  mutate(DBH = D4/10) %>%
  dplyr::select(-D4) %>%
  filter(!is.na(census),
         !is.na(patch))

censuses <-
  data.selected %>% mutate(
    DBH_group = case_when(
      DBH < 10 ~ 0,
      DBH < 20 ~ 1,
      DBH < 30 ~ 2,
      DBH < 40 ~ 3,
      DBH < 50 ~ 4,
      DBH < 60 ~ 5,
      DBH < 70 ~ 6,
      DBH < 80 ~ 7,
      DBH < 90 ~ 8,
      DBH < 100 ~ 9,
      TRUE ~ 10)) %>%
      mutate(PFT = classify.sp(species = species) %>% pull(PFT))

census.sum <- censuses %>%
  group_by(census,PFT,patch,DBH_group) %>%
  summarise(N = n(),
            .groups = "keep") %>%
  ungroup() %>%
  complete(census = 1:4,
           patch = 1:25,
           PFT = c(2,3,4),
           DBH_group = seq(0,10),
           fill = list(N = 0)) %>%
  group_by(census,PFT,DBH_group) %>%
  summarise(sd = sd(N, na.rm = TRUE),
            N.m = mean(N,na.rm = TRUE),
            N = length(DBH_group),
            .groups = "keep")

tot.sum <- censuses %>%
  group_by(census,patch,DBH_group) %>%
  summarise(N = n(),
            .groups = "keep") %>%
  ungroup() %>%
  complete(census = 1:4,patch = 1:25,DBH_group = seq(0,10),fill = list(N = 0)) %>%
  group_by(census,DBH_group) %>%
  summarise(sd = sd(N, na.rm = TRUE),
            N.m = mean(N,na.rm = TRUE),
            N = length(DBH_group),
            se = sd/sqrt(N),
            .groups = "keep")

fac = 10000/20/20

ggplot() +
  geom_errorbar(data = tot.sum %>%
                  filter(DBH_group > 0),
                mapping = aes(x = as.factor(DBH_group),
                              ymin = 0.75*(N.m)*fac, ymax = (N.m + sd)*fac),
                width = 0.2) +
  geom_bar(data = census.sum %>%
             filter(DBH_group > 0),
           mapping = aes(x = as.factor(DBH_group),
                         y = N.m*fac, fill = as.factor(PFT)),
           stat="identity",
           position="stack") +
  labs(fill = "Plant Functional Type", x = "DBH [cm]", y = "Tree density [ind/ha]") +
  scale_fill_manual(values = c("#9FFF8C","#44CC29","#137300"),
                    labels = c("Early","Mid","Late")) +
  scale_x_discrete(breaks = 1:10,
                     labels = c("10-20","20-30","30-40","40-50",
                                "50-60","60-70","70-80","80-90",
                                "90-100", ">100") ) +
  theme_bw() +
  facet_wrap(~census) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = c(0.85, 0.85),
        text = element_text(size = 24))


################################################################################


dir.patch <- file.path(getwd(),
                       "/inputs/csspss/")

census.target <- 1

for (iplot in seq(1,nrow(data.plot))){

  cplotID <- data.plot$PlotID[iplot]
  cplot.code <- data.plot$`Plot Code`[iplot]
  cdir <- file.path(dir.patch,cplot.code)

  dir.create(cdir,showWarnings = FALSE)

  clat <- data.census.interval %>%
    filter(PlotCode == cplot.code) %>%
    pull(Lat) ;

  clon <- data.census.interval %>%
    filter(PlotCode == cplot.code) %>%
    pull(Lon)

  csspssfile_name <- paste0(cplot.code,
                            ".lat",sprintf("%.3f",clat),
                            ".lon",sprintf("%.3f",clon))

  plot.patch <- censuses %>%
    filter(census == census.target,
           PlotID == cplotID) %>%
    group_by(patch) %>%
    summarise(trk = 2,
              time = data.census %>%
                filter(census.no == census.target) %>%
                pull(census.date),
              age = 0,water = 0,fsc = 0.15,
              stsc = 6,stsl = 6,ssc = 4.5,
              lai = 5,msn = 0.6,fsn = 0.003,
              nep = 0,gpp = 0,rh = 0,
              .groups = "keep") %>%
    arrange(patch) %>%
    ungroup() %>%
    mutate(area = 1/length(unique(patch))) %>%
    dplyr::select(time,patch,trk,age,area,water,fsc,stsc,stsl,ssc,lai,msn,fsn,nep,gpp,rh)

  write_pss(pss = plot.patch,path_prefix = file.path(dir.patch,cplot.code,cplot.code),latitude = clat,longitude = clon)

  plot.ind <- censuses %>%
    filter(census == census.target,
           PlotID == cplotID) %>%
    arrange(patch) %>%
    rename(pft = PFT,
           dbh = DBH) %>%
    ungroup() %>%
    mutate(hite = 0,
           bdead = 0,
           balive = 0,
           lai = 0,
           n = 1/(20*20),
           time = data.census %>%
             filter(census.no == census.target) %>%
             pull(census.date),
           cohort = 1:length(TreeID)) %>%
    dplyr::select(time,patch,cohort,dbh,hite,pft,n,bdead,balive,lai)


  write_css(css = plot.ind,path_prefix = file.path(dir.patch,cplot.code,cplot.code),latitude = clat,longitude = clon)

  system2("rsync",
          c("-avz",
            file.path(dir.patch,cplot.code),
            "hpc:/data/gent/vo/000/gvo00074/ED_common_data/inits/"))
}


