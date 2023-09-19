###created by Huiyun Wu
###modified by Huiyun Wu 9/14/23

# install.packages(c("ggplot2", "ggpubr", "tidyverse", "broom", "AICcmodavg","NADA","visdat","naniar"))



library(tidyverse)
library(tibble)
library(dplyr)
library(ggplot2)
library(dplyr)

setwd("/Users/huiyunwu/Desktop/Virus_particle/")
#load stacked dataframe from row data#  
# df_stacked0<-read.csv(file="processed_data/PAV_complete_stacked0.csv")
# summary(df_stacked0)

s3tos32<-read_csv(file = "data/processed_data/S3to32_back.cal_dup.csv")
colnames(s3tos32)
s3tos32

s3tos32 %>% 
  arrange(run)

df<-s3tos32 %>% 
  arrange(run) %>% 
  mutate(run = tolower(run),
         run = str_replace(run, "adenovirus","adenovirus"),
         target = str_sub(run, start =1, end =3),
         target = str_replace(target, "cra","crAssphage56"),
         target = str_replace(target, "ent", "enterovirus"),
         target = str_replace(target, "ade", "adenovirus"),
         target = str_replace(target, "pmm", "pmmov"),
         target = str_replace(target, "crAssphage56", "crAssphage56"),
         target = str_replace(target, "nor", "norovirus")) 



df %>% 
  filter(qualitycheck == "Pass") %>% 
  group_by(name,target,wwtp, collection.date,pore.size) %>%
  summarize(mean =mean(finalcon.cp.L, na.rm = T)) %>%
  drop_na(wwtp) %>% 
  arrange(desc(target)) %>% 
  slice_tail(n=10)

df2 %>% 
  arrange(wwtp) %>% 
  slice_max(wwtp, n=3)

df2 %>% 
  arrange(desc(wwtp))

df2 %>% 
  tail()
  
df2 %>% 
  head()

long <-df %>% 
  filter(qualitycheck == "Pass") %>% 
  group_by(name,target,wwtp, vol.L, collection.date,pore.size) %>%
  summarize(cp.per.L =mean(finalcon.cp.L, na.rm = T)) %>%
  drop_na(wwtp) %>% 
  # arrange(collection.date) %>% 
  arrange(target)

Vol.<-read_csv(file = "data/processed_data/Vol.Sample.csv",
               col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                                name = col_character()))

wide<-pivot_wider(long, id_cols = c(name,wwtp, collection.date, pore.size, vol.L), names_from = target, values_from = cp.per.L)
# write_csv(wide, file = "data/processed_data/wide.csv")


# unstacked<-left_join(Vol.[,c(1:8)],wide) %>% 
#   arrange(crAssphage56)
# write_csv(unstacked, file = "data/processed_data/unstacked.csv")
  
 # write_csv(wide, file = "data/processed_data/wide.csv")
unstacked<- read_csv(file = "data/processed_data/unstacked_matched_colored_v2.csv",
                     col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                                      name = col_character()))
summary(unstacked)

unstacked %>% 
  group_by(wrf) %>% 
  summarise(count = n())


#####count for pmmov, table 3#####
unstacked %>% 
  filter(pmmov > 15) %>% 
  group_by(extn.type) %>% 
  summarise(count = n())


unstacked %>% 
  filter(pmmov > 15) %>% 
  summarise(count = n())/137

unstacked %>% 
  filter(crAssphage56 > 15) %>% 
  summarise(count = n())/148

unstacked %>% 
  filter(adenovirus > 15) %>% 
  summarise(count = n())/148

unstacked %>% 
  filter(norovirus > 15) %>% 
  summarise(count = n())/148

unstacked %>% 
  filter(enterovirus > 15) %>% 
  summarise(count = n())/148

unstacked %>% 
  filter(ms2 > 15) %>% 
  summarise(count = n())/94

unstacked %>% 
  filter(crAssphage56 > 15) %>% 
  summarise(count = n())/148


#count for ms2, table 3
unstacked[c(1:94),] %>% 
  group_by(extn.type) %>% 
  summarise(count = n())
  
unstacked %>% 
  group_by(extn.type) %>% 
  summarise(count = n()) 


#####supplimentary table s3#####
filter_100 <- unstacked %>% 
  select(pore.size, crAssphage56, pmmov, adenovirus,norovirus,enterovirus) %>% 
  filter(pore.size == "100 um") 

 summary(filter_100)

filter_20 <- unstacked %>% 
   select(name, pore.size, crAssphage56, pmmov, adenovirus,norovirus,enterovirus) %>% 
   filter(pore.size == "20 um") 
 
 summary(filter_20)
 

filter_3 <- unstacked %>% 
   select(name, pore.size, crAssphage56, pmmov, adenovirus,norovirus,enterovirus) %>% 
   filter(pore.size == "3 um") 
 

 summary(filter_3)
 
 
filter_0.45 <- unstacked %>% 
   select(pore.size, crAssphage56, pmmov, adenovirus,norovirus,enterovirus) %>% 
   filter(pore.size == "0.45 um") 
 
print(summary(filter_0.45))
 

filtrate <- unstacked %>% 
   select(pore.size, crAssphage56, pmmov, adenovirus,norovirus,enterovirus) %>% 
   filter(pore.size == "0.45 um filtrate") 
 
 summary(filtrate)
 
#######load unstacked dataframe #####
colnames(unstacked)
summary(unstacked)
# create cleaner dataframe from unstacked dataframe
# stacked<-unstacked %>%
#   pivot_longer(
#     cols = crAssphage56:pmmov,
#     names_to = "target",
#     values_to = "finalcon.cp.L")
# write.csv(stacked, file="/Users/huiyunwu/Desktop/Virus_particle/data/processed_data/PAV_stacked_v2.csv",row.names = F)

#load stacked dataframe
stacked<-read_csv(file="data/processed_data/pav_stacked_v2.csv")
colnames(stacked)
###### load long dataframe############
stacked$finalcon.cp.L<-stacked$finalcon.cp.L+1
stacked$log_targets<-log10(stacked$finalcon.cp.L)

###colors###
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(RColorBrewer)
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 6, name = 'Accent')


unstacked<-read_csv(file="data/processed_data/unstacked_matched_colored_v2.csv")
data1<-read_csv(file="data/processed_data/unstacked_matched_colored_v2.csv")
colnames(data1)
data1$ms2<-data1$ms2+1
data1$adenovirus<-data1$adenovirus+1
data1$enterovirus<-data1$enterovirus+1
data1$norovirus<-data1$norovirus+1
data1$crAssphage56<-data1$crAssphage56+1
data1$pmmov<-data1$pmmov+1
log_targets<-log10(data1[,c("crAssphage56","ms2","pmmov","adenovirus","enterovirus","norovirus")])
colnames(log_targets)
env <- data1[,c("name","wwtp","wrf", "collection.date","pore.size","extn.type","season","vol.L")]
log_data2<-cbind(env, log_targets)
str(log_data2)
# log_data2$name<-as.character(log_data2$name)
head(log_data2)

##rank the factors##
#Pore.size
stacked$pore.size<-factor(stacked$pore.size, levels=c('100 um','20 um','3 um','0.45 um','0.45 um filtrate'),ordered = TRUE)
#Target
stacked$target<-factor(stacked$target,
                        levels=c('crAssphage56','pmmov','ms2','adenovirus','norovirus','enterovirus'),ordered = TRUE)
#season
stacked$season<-factor(stacked$season,
                        levels=c('spring','summer','fall','winter'),ordered = TRUE)
colnames(stacked)

######by pore size#######
ggplot(stacked,
       aes(x = pore.size,
           y=log_targets,
           fill=target)) +
  # geom_boxplot()+
   geom_violin()+
  geom_boxplot(width=0.1) +
  facet_wrap(. ~target, ncol = 2, nrow =3)+
  geom_point(position = position_jitter(seed = 1, width = 0.2), size =0.3) +
  scale_fill_brewer(palette = "Accent")+
  # scale_fill_manual(values=cbPalette)+
  labs(title="Viruses in secondary effluent",
       y= "Target Concentration (logGC/L)",
       x="Membrane pore size")+
  theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        # axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(color = "Black",face = "bold", size = rel(1.5)),)
# ggsave("plots/By_poresize_2.pdf")

###POST-HOC TEST COMPARISON Pore sizes##
unstacked<-read_csv(file="data/processed_data/unstacked_matched_colored_v2.csv")
data1<-read_csv(file="data/processed_data/unstacked_matched_colored_v2.csv")
colnames(data1)
data1$ms2<-data1$ms2+1
data1$adenovirus<-data1$adenovirus+1
data1$enterovirus<-data1$enterovirus+1
data1$norovirus<-data1$norovirus+1
data1$crAssphage56<-data1$crAssphage56+1
data1$pmmov<-data1$pmmov+1
log_targets<-log10(data1[,c("crAssphage56","ms2","pmmov","adenovirus","enterovirus","norovirus")])
colnames(log_targets)
env <- data1[,c("name","wwtp","wrf", "collection.date","pore.size","extn.type","season","vol.L")]
log_data2<-cbind(env, log_targets)
str(log_data2)
# log_data2$name<-as.character(log_data2$name)
head(log_data2)

crAssphage56_pore.size<-aov(log_data2$crAssphage56~log_data2$pore.size)
summary(crAssphage56_pore.size)
TukeyHSD(crAssphage56_pore.size)

pmmov_pore.size<-aov(log_data2$pmmov~log_data2$pore.size)
summary(pmmov_pore.size)
TukeyHSD(pmmov_pore.size)

adenovirus_pore.size<-aov(log_data2$adenovirus~log_data2$pore.size)
summary(adenovirus_pore.size)
TukeyHSD(adenovirus_pore.size)  

enterovirus_pore.size<-aov(log_data2$enterovirus~log_data2$pore.size)
summary(enterovirus_pore.size)
TukeyHSD(enterovirus_pore.size)

norovirus_pore.size<-aov(log_data2$norovirus~log_data2$pore.size)
summary(norovirus_pore.size)
TukeyHSD(norovirus_pore.size)

######By extraction type#######
ggplot(stacked,
       aes(x = extn.type,
           y=log_targets,
           fill=target)) +
  # geom_boxplot()+
  geom_violin()+
  facet_wrap(. ~target, ncol = 2, nrow = 3)+
  geom_point(position = position_jitter(seed = 1, width = 0.3), size =0.3) +
  scale_fill_brewer(palette = "Accent")+
  # scale_fill_manual(values=cbPalette)+
  labs(title=NULL,
       y= "Target Concentration (logGC/L)",
       x="Extraction type")+
  theme(strip.background = element_rect(fill = "lightgrey"),
        strip.text = element_text(color = "Black",face = "bold", size = rel(1.5)),)
# ggsave("plots/By_EnxtType2.pdf")


###POST-HOC TEST COMPARISON extn.type###
crAssphage56_extn.type<-aov(log_data2$crAssphage56~log_data2$extn.type)
summary(crAssphage56_extn.type)
TukeyHSD(crAssphage56_extn.type)


adenovirus_extn.type<-aov(log_data2$adenovirus~log_data2$extn.type)
summary(adenovirus_extn.type)
TukeyHSD(adenovirus_extn.type)  

enterovirus_extn.type<-aov(log_data2$enterovirus~log_data2$extn.type)
summary(enterovirus_extn.type)
TukeyHSD(enterovirus_extn.type)

norovirus_extn.type<-aov(log_data2$norovirus~log_data2$extn.type)
summary(norovirus_extn.type)
TukeyHSD(norovirus_extn.type)

pmmov_extn.type<-aov(log_data2$pmmov~log_data2$extn.type)
summary(pmmov_extn.type)
TukeyHSD(pmmov_extn.type)






######by WWTP######
ggplot(stacked,
       aes(x = wrf,
           y=log_targets,
           fill=target)) +
  geom_boxplot()+
  # geom_violin()+
  # facet_wrap(. ~target, ncol = 2, nrow = 3)+
  geom_point(position = position_jitter(seed = 1, width = 0.2), size =0.3) +
  scale_fill_brewer(palette = "Accent")+
  # scale_fill_manual(values=cbPalette)+
  labs(title=NULL,
       y= "Target Concentration (logGC/L)",
       x="WWTP")+
  theme_bw()
 # ggsave("plots/By_WWTP_1.pdf")
#POST-HOC TEST COMPARISON wrf##
crAssphage56_wrf<-aov(log_data2$crAssphage56~log_data2$wrf)
summary(crAssphage56_wrf)
TukeyHSD(crAssphage56_wrf)

pmmov_wrf<-aov(log_data2$pmmov~log_data2$wrf)
summary(pmmov_wrf)
TukeyHSD(pmmov_wrf)

adenovirus_wrf<-aov(log_data2$adenovirus~log_data2$wrf)
summary(adenovirus_wrf)
TukeyHSD(adenovirus_wrf)  

enterovirus_wrf<-aov(log_data2$enterovirus~log_data2$wrf)
summary(enterovirus_wrf)
TukeyHSD(enterovirus_wrf)

norovirus_wrf<-aov(log_data2$norovirus~log_data2$wrf)
summary(norovirus_wrf)
TukeyHSD(norovirus_wrf)


########By season#######
ggplot(stacked,
       aes(x = season,
           y=log_targets,
           fill=target)) +
  geom_boxplot(width=1)+
  geom_point(position = position_jitter(seed = 1, width = 0.3), size =0.3) +
  # facet_wrap(. ~season, ncol = 2, nrow = 3)+
  scale_fill_brewer(palette = "Accent")+
  # scale_fill_manual(values=cbPalette)+
  labs(title= NULL,
       y= "Target Concentration (logGC/L)",
       x="season")+
  theme_bw()
 # ggsave("plots/By_season1.pdf")

##POST-HOC TEST COMPARISON season##
crAssphage56_season<-aov(log_data2$crAssphage56~log_data2$season)
summary(crAssphage56_season)
TukeyHSD(crAssphage56_season)

pmmov_season<-aov(log_data2$pmmov~log_data2$season)
summary(pmmov_season)
TukeyHSD(pmmov_season)

adenovirus_season<-aov(log_data2$adenovirus~log_data2$season)
summary(adenovirus_season)
TukeyHSD(adenovirus_season)  

enterovirus_season<-aov(log_data2$enterovirus~log_data2$season)
summary(enterovirus_season)
TukeyHSD(enterovirus_season)

norovirus_season<-aov(log_data2$norovirus~log_data2$season)
summary(norovirus_season)
TukeyHSD(norovirus_season)






#####unstacked data##
# data1<-read.csv(file="processed_data/filters_47mm.data_unstacked_v2.csv")
# 
# colnames(data1)
# 
# data1$adenovirus<-as.numeric(data1$adenovirus)
# data1$enterovirus<-as.numeric(data1$enterovirus)
# data1$norovirus<-as.numeric(data1$norovirus)
# data1$crAssphage56<-as.numeric(data1$crAssphage56)
# data1$ms2<-as.numeric(data1$ms2)
# 
# data1$adenovirus[is.na(data1$adenovirus)] <- 0
# data1$enterovirus[is.na(data1$enterovirus)] <- 0
# data1$norovirus[is.na(data1$norovirus)] <- 0
# data1$crAssphage56[is.na(data1$crAssphage56)] <- 0
# data1$ms2[is.na(data1$ms2)] <- 0

stacked<-read_csv(file="data/processed_data/PAV_stacked_v2.csv",
               col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                                name = col_character()))


colnames(stacked)
######check if there is duplicate in data##
# df_stacked %>%
#   dplyr::group_by(Name, WWTP, wrf, Collection.date, Pore.size, Vol.L, Target) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L)

# df_wide <- pivot_wider(df_stacked, id_cols = c(Name,WWTP, wrf, Collection.date, Pore.size, extn.type,season, Vol.L), names_from = Target, values_from = finalcon.cp.L)
# colnames(df_wide)
# 
#  write_csv(df_wide,file="processed_data/PAV_complete_unstacked.csv")


######load wide dataframe######
unstacked<-read_csv(file="data/processed_data/unstacked_matched_colored_v2.csv",
                  col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                                   name = col_character()))

unstacked
######reformat the data for log10 #####
data1<-unstacked
data1$ms2<-data1$ms2+1
data1$adenovirus<-data1$adenovirus+1
data1$enterovirus<-data1$enterovirus+1
data1$norovirus<-data1$norovirus+1
data1$crAssphage56<-data1$crAssphage56+1
data1$pmmov<-data1$pmmov+1
######prepare log10 data frame############
colnames(data1)
log_targets<-log10(data1[,c("crAssphage56","ms2","pmmov","adenovirus","enterovirus","norovirus")])
colnames(log_targets)
env <- data1[,c("name","wwtp","wrf", "collection.date","pore.size","extn.type","season","vol.L")]
log_data2<-cbind(env, log_targets)
str(log_data2)
# log_data2$name<-as.character(log_data2$name)
head(log_data2)

###qq plot check normal distribution######
# qqnorm(log_data2$crAssphage56)
# qqline(log_data2$crAssphage56)
library("car")
qqPlot(log_data2$crAssphage56)
qqPlot(log_data2$pmmov)
qqPlot(log_data2$adenovirus)
qqPlot(log_data2$enterovirus)
qqPlot(log_data2$norovirus)



#####one way ANOVA #####
#on pore size##
# fit_crAssphage <- lm(crAssphage56 ~ pore.size, data=log_data2)
# summary(fit_crAssphage)
# fit_pmmov<-lm(pmmov ~ pore.size, data=log_data2)
# summary(fit_pmmov)
# fit_adenovirus <- lm(adenovirus ~ pore.size, data=log_data2)
# summary(fit_adenovirus)
# fit_enterovirus <- lm(enterovirus ~ pore.size, data=log_data2)
# summary(fit_enterovirus)
# fit_norovirus <- lm(norovirus ~ pore.size, data=log_data2)
# summary(fit_norovirus)

#one way ANOVA on wrf#
# fit_crAssphage56 <- lm(crAssphage56 ~ wrf, data=log_data2)
# summary(fit_crAssphage56)
# fit_pmmov<-lm(pmmov ~ wrf, data=log_data2)
# summary(fit_pmmov)
# fit_adenovirus <- lm(adenovirus ~ wrf, data=log_data2)
# summary(fit_adenovirus)
# fit_enterovirus <- lm(enterovirus ~ wrf, data=log_data2)
# summary(fit_enterovirus)
# fit_norovirus <- lm(norovirus ~ wrf, data=log_data2)
# summary(fit_norovirus)

#one way ANOVA on season#
# fit_crAssphage56 <- lm(crAssphage56 ~ season, data=log_data2)
# summary(fit_crAssphage56)
# fit_pmmov<-lm(pmmov ~ season, data=log_data2)
# summary(fit_pmmov)
# fit_adenovirus <- lm(adenovirus ~ season, data=log_data2)
# summary(fit_adenovirus)
# fit_enterovirus <- lm(enterovirus ~ season, data=log_data2)
# summary(fit_enterovirus)
# fit_norovirus <- lm(norovirus ~ season, data=log_data2)
# summary(fit_norovirus)
# 


#one way ANOVA on extraction type#
# fit_crAssphage56 <- lm(crAssphage56 ~ extn.type, data=log_data2)
# summary(fit_crAssphage56)
# fit_ms2 <- lm(ms2 ~ extn.type, data=log_data2)
# summary(fit_ms2)
# fit_adenovirus <- lm(adenovirus ~ extn.type, data=log_data2)
# summary(fit_adenovirus)
# fit_enterovirus <- lm(enterovirus ~ extn.type, data=log_data2)
# summary(fit_enterovirus)
# fit_norovirus <- lm(norovirus ~ extn.type, data=log_data2)
# summary(fit_norovirus)
# fit_pmmov<-lm(pmmov ~ extn.type, data=log_data2)
# summary(fit_pmmov)


#
#Winter 12-2
#Spring 3-5
#Summer 6-8
#Fall 9-11

###

####correlation between two targets####
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}
ggplotRegression(lm(crAssphage56~adenovirus, data=unstacked))

######multiple testing adjustment p value, spearman test ########
library(Hmisc)
##set detection limit to 15 cp/L
unstacked[is.na(unstacked)]=0.59 ###replace NA to 0.59 log10/L, half of detection limit 15
matrix1<-rcorr(as.matrix(unstacked[,c("crAssphage56","pmmov","adenovirus","norovirus","enterovirus")]), type = "spearman")
print(matrix1)
 capture.output(print(matrix1), file = "data/processed_data/spearman_corr2.csv")
# write.table(as.data.frame(matrix1),file="spearman_corr.csv", quote=F,sep=",",row.names=F)
 
#####non-metric multidimensional scaling NMDS ######
data3<- read_csv(file = "data/processed_data/unstacked_matched_colored_v2.csv",
                      col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                                       name = col_character()))
 
data3[data3 == 0] <- NA
#remove ms2
colnames(data3)
data3<-data3 %>% 
  select(!ms2)
data3[is.na(data3)]=0.59 ###replace NA to 0.65 log10/L, half of detection limit 32.4
library(readr)
library(vegan)
m_data3<-data3
# m_data3=as.matrix(data3)
 set.seed(123)
 # nmds = metaMDS(m_data3, distance = "bray")
# 'comm' has negative data: 'autotransform', 'noshare' and 'wascores' set to FALSE
# Error in distfun(comm, method = distance, ...) :
  # input data must be numeric
head(data3)
colnames(data3)
virus = data3[,c(9:13)]
m_virus = as.matrix(virus)
m_virus<-data3 %>% 
  select(crAssphage56, pmmov, adenovirus, enterovirus, norovirus)
set.seed(123)
nmds = metaMDS(m_virus, distance = "bray")
nmds = metaMDS(m_virus, distance = "bray")
nmds
# Call:
#   metaMDS(comm = m_virus, distance = "bray") 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(sqrt(m_virus)) 
# Distance: bray 
# 
# Dimensions: 2 
# Stress:     0.1376374 
# Stress type 1, weak ties
# Best solution was repeated 4 times in 20 tries
# The best solution was from try 5 (random start)
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(m_virus))’ 


nmds$stress

test<-scores(nmds)

# scores(nmds) %>%
#   as_tibble(rownames = "samples") %>%
#   inner_join(., sample_lookup, by="samples") %>%
#   mutate(period = if_else(day < 10, 'early', 'late')) %>%
#   ggplot(aes(x=NMDS1, y=NMDS2, color=period)) +
#   geom_point()



plot(nmds)
as.data.frame(scores(nmds)$site)
data.scores = as.data.frame(scores(nmds)$site)
# View(data3)
data.scores$wrf = data3$wrf
data.scores$date = data3$Collection.date
data.scores$pore.size = data3$pore.size
data.scores$season = data3$season
data.scores$extn.type = data3$extn.type
# data.scores$Influent = pathogen$Influent
head(data.scores)
library(ggplot2) 
#by pore size
# xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
#  +   geom_point(size = 4, aes( shape = pore.size, colour = Pore. size))+
#  +   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
#            +         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
#            +         legend.text = element_text(size = 12, face ="bold", colour ="black"),
#            +         legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
#            +         axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
#            +         legend.title = element_text(size = 14, colour = "black", face = "bold"),
#            +         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
#            +         legend.key=element_blank()) +
#  +   labs(x = "NMDS1", colour = "Time", y = "NMDS2", shape = "Type")  +
#  +   scale_colour_manual(values = c("#009E73", "#E69F00"))
data.scores$pore.size<-factor(data.scores$pore.size, levels=c('100 um','20 um','3 um','0.45 um','0.45 um filtrate'),ordered = TRUE)

xx = ggplot(data.scores, 
            aes(x = NMDS1, y = NMDS2, color= pore.size, fill=pore.size)) + 
  stat_ellipse(geom="polygon",level=0.9, alpha=0.3,show.legend = FALSE)+
  geom_point(size = 1)+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  scale_colour_manual(values = c("#009E73", "#E69F00","#FF0000","gray45", "#2166AC")
  )+
  scale_fill_manual(values = c("palegreen", "wheat","lightpink","gray85","#4393C3")
  ) 

xx
#####R markdown code####
library(rmarkdown)
library(tinytex)
render("PAV_writing_rmd.Rmd",output_format = "all")
