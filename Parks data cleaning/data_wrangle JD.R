
# NOTES:
# (fixed some column names in csv)
# should fix there being two "Plot" columns in data (generally bad practice)
# changed WIR column name to match others TOTAL FR
# took out "-" from height column in data
# added plot and subplot col to alones
# changed alone stemcount to stems
# added a "park" column to neighborhood sheets;  made "plot" and subplot lowercase
# erased "skip" from WIR neighborhood file under PLECON column.   no other data for those neigh, and was causing that column to be read as "chr"

library(lme4)
library(tidyverse)
library(stringr)
theme_set(theme_minimal())

# loading in data (both focal and neighborhood data)----
rf_focals <- read.csv("data_raw_JD/riverfront focals.csv", na.strings=c("NA","NaN", " ", "") )
rf_neighborhoods <- read.csv("data_raw_JD/riverfront neighborhoods.csv", na.strings=c("NA","NaN", " ", "") )

wir_focals <- read.csv("data_raw_JD/wir focals.csv", na.strings=c("NA","NaN", " ", "") )
wir_neighborhoods <- read.csv("data_raw_JD/wir neighborhoods.csv", na.strings=c("NA","NaN", " ", "") )

br_focals <- read.csv("data_raw_JD/br focals.csv", na.strings=c("NA","NaN", " ", "") )
br_neighborhoods <- read.csv("data_raw_JD/br neighborhoods.csv", na.strings=c("NA","NaN", " ", "") )

sem_focals <- read.csv("data_raw_JD/sem focals.csv", na.strings=c("NA","NaN", " ", "") )
sem_neighborhoods <- read.csv("data_raw_JD/sem neighborhoods.csv", na.strings=c("NA","NaN", " ", ""))

head(sem_focals,2); dim(sem_focals)
head(sem_neighborhoods,2); dim(sem_neighborhoods)

table(sem_neighborhoods$focal_collect)

head(sem_focals,2)
head(wir_focals,2)
head(br_focals,2)

alones_all <- read.csv("data_raw_JD/alones.csv")
names(alones_all) <- tolower(names(alones_all))
alones_all <- alones_all %>% filter(!grepl('GONE', notes),
                                    !grepl('DEAD', notes),
                            !grepl("couldn't find", notes))

head(alones_all); dim(alones_all)


# Put together Non-riverfront parks 
# dropping all rows with notes Gone
focals <- rbind(sem_focals, br_focals, wir_focals)
focals <- focals %>% filter(!grepl('GONE', NOTES),
                                    !grepl("couldn't find", NOTES))#,
                                    # !grepl('BT', NOTES)) 
names(focals) <- tolower(names(focals))
head(focals,2); dim(focals)
str(focals)

# COLLOM ----

##  Eugene Parks ----
# Fitness = # seeds = [ (# top infl x # flrs/big) + (# side infl x # flrs/small) ]  * 3 seeds per fruit 
# TOTALFR is the total number of infloresences 
# TOPFR is number of fruits on the top infloresence
# SIDEFR is the number of fruits on side infl
collom_parks <- filter(focals, species=="COLLOM") %>% rowwise() %>%
  mutate(stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr), # was top_infl=stems but changed because there were some with stems > total.fr (stems without infl)
         fruits=topfr*top_infl + sidefr*side_infl)
print(arrange(collom_parks,desc(fruits)), n=8, width=Inf); dim(collom_parks)

## Riverfront ----
# TOTALFR is the total number of infloresences 
# TOPFR is number of fruits on the top infloresence
# SIDEFR is the number of fruits on side infl
# Number of side infl = fr1 – stemcount 

rf_focals <- rf_focals %>% filter(!grepl('GONE', NOTES),
                                  !grepl('DEAD', NOTES),
                                  !grepl("couldn't find", NOTES))#,
names(rf_focals) <- tolower(names(rf_focals))

head(filter(rf_focals, species=="COLLOM"))
collom_rf <- filter(rf_focals, species=="COLLOM")
collom_rf <- collom_rf %>% rowwise() %>%
  mutate(total.fr = rowSums(across(fr1:fr5), na.rm=TRUE),
         stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr),
         fruits=topfr*top_infl + sidefr*side_infl) %>% 
  ungroup()
print(arrange(collom_rf, desc(fruits)), n=8, width=Inf)

# so collom_rf and collom_parks have fruit estimates for plants with fruit counts

## get fruit~height relationships ----
# need to get estimates for the plants without fruit counts based on height  
# (doing this separately for 'diverse' plot plants versus alone plants)

# Put together diverse Collom (parks + rf)
collom_all <- rbind(dplyr::select(collom_rf, park:stems, biomass:notes, fruits), dplyr::select(collom_parks, park:stems, biomass:notes, fruits)) %>%
  as.data.frame()
head(collom_all, n=8, width=Inf)

ggplot(collom_all, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="COLLOM")
ggsave("Figures/collom_height_fruits_park.pdf", width=6, height=6)

ggplot(collom_all, aes(x=height, y=fruits, color=plot))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="COLLOM")
ggsave("Figures/collom_height_fruits.pdf", width=6, height=6)

ggplot(collom_all, aes(x=height, y=biomass, color=plot))+
  geom_point() + geom_smooth(method='lm')+
  labs(title="COLLOM")
ggsave("Figures/collom_height_biomass.pdf", width=6, height=6)


collom_all_sub <- collom_all %>% filter(!(is.na(fruits)))  # is.na(height) &
head(collom_all_sub, n=4, width=Inf); dim(collom_all_sub)

# fitting park-specific relationships
m1 <- glm(fruits ~ height*park, data=collom_all_sub, family=poisson(link="log"))
summary(m1)

collom_all$fruits_pred <- ceiling(predict(m1, collom_all, type = "response"))
head(collom_all, n=8, width=Inf)
collom_all
ggplot(collom_all, aes(x=fruits, y=fruits_pred))+
  geom_point() + geom_smooth(method='lm')+
  geom_abline(intercept=0,slope=1, linetype='dashed')

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
collom_all <- collom_all %>%   # don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))


## alone Collom ----

collom_alone <- filter(alones_all, species=="COLLOM") %>% 
  rowwise() %>%
  mutate(total.fr = rowSums(across(fr1:fr5), na.rm=TRUE),
         stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr),
         fruits=topfr*top_infl + sidefr*side_infl)  %>%
  ungroup()%>% 
  as.data.frame()
  
head(collom_alone, n=3, width=Inf)

ggplot(collom_alone, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="COLLOM alones")
ggsave("Figures/collom_ALONES_height_fruits_byPark.pdf", width=6, height=6)

ggplot(collom_alone, aes(x=height, y=fruits, color=park))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="COLLOM alones")
ggsave("Figures/collom_ALONES_height_fruits.pdf", width=6, height=6)



collom_alone_sub <- collom_alone %>% filter(!(is.na(fruits)))  # fit fruit ~ height relationship with those plants with both
head(collom_alone_sub, n=4, width=Inf); dim(collom_all_sub)

# do these regressions by park
m2 <- glm(fruits ~ height*park , data=collom_alone_sub, family=poisson(link="log"))
summary(m2)

collom_alone$fruits_pred <- predict(m2, collom_alone, type = "response")   # type = "response"    needed to get fruits, not log(fruits)

ggplot(collom_alone, aes(x=fruits, y=fruits_pred))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
head(collom_alone); dim(collom_alone)
collom_alone <- collom_alone %>%   # checked and don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))%>% 
  dplyr::select(park:height, biomass:notes, fruits) %>% 
  as.data.frame()
collom_all <- collom_all %>% select(-c(stems,fruits_pred))
  
head(collom_alone, n=4, width=Inf); dim(collom_alone)
head(collom_all, n=4, width=Inf); dim(collom_all)


# NAVSQU ----

##  Eugene Parks ----
# Fitness = # seeds = [ (# top infl x # flrs/big) + (# side infl x # flrs/small) ]  * 3 seeds per fruit 
# TOTALFR is the total number of infloresences 
# TOPFR is number of fruits on the top infloresence
# SIDEFR is the number of fruits on side infl
navsqu_parks <- filter(focals, species=="NAVSQU") %>% rowwise() %>%
  mutate(stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr), # was top_infl=stems but changed because there were some with stems > total.fr (stems without infl)
         fruits=topfr*top_infl + sidefr*side_infl) %>%
  ungroup()
print(arrange(navsqu_parks,desc(fruits)), n=8, width=Inf); dim(navsqu_parks)

## Riverfront ----
# TOTALFR is the total number of infloresences 
# TOPFR is number of fruits on the top infloresence
# SIDEFR is the number of fruits on side infl
# Number of side infl = fr1 – stemcount 

rf_focals <- rf_focals %>% filter(!grepl('GONE', notes),
                                  !grepl('DEAD', notes),
                                  !grepl("couldn't find", notes))#,
names(rf_focals) <- tolower(names(rf_focals))

head(filter(rf_focals, species=="NAVSQU"))
navsqu_rf <- filter(rf_focals, species=="NAVSQU")
navsqu_rf <- navsqu_rf %>% rowwise() %>%
  mutate(total.fr = rowSums(across(fr1:fr5), na.rm=TRUE),
         stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr),
         fruits=topfr*top_infl + sidefr*side_infl) %>%
  ungroup()
print(arrange(navsqu_rf, desc(fruits)), n=8, width=Inf)

# so navsqu_rf and navsqu_parks have fruit estimates for plants with fruit counts

## get fruit~height relationships ----
# need to get estimates for the plants without fruit counts based on height  
# (doing this separately for 'diverse' plot plants versus alone plants)

# Put together diverse navsqu (parks + rf)
navsqu_all <- rbind(dplyr::select(navsqu_rf, park:stems, biomass:notes, fruits), dplyr::select(navsqu_parks, park:stems, biomass:notes, fruits)) %>%
  as.data.frame()
head(navsqu_all, n=8, width=Inf)

ggplot(navsqu_all, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="navsqu")
ggsave("Figures/navsqu_height_fruits_park.pdf", width=6, height=6)

ggplot(navsqu_all, aes(x=height, y=fruits, color=plot))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="navsqu")
ggsave("Figures/navsqu_height_fruits.pdf", width=6, height=6)

ggplot(navsqu_all, aes(x=height, y=biomass, color=plot))+
  geom_point() + geom_smooth(method='lm')+
  labs(title="navsqu")
ggsave("Figures/navsqu_height_biomass.pdf", width=6, height=6)


navsqu_all_sub <- navsqu_all %>% filter(!(is.na(fruits)))  # is.na(height) &
head(navsqu_all_sub, n=4, width=Inf); dim(navsqu_all_sub)

m1 <- glm(fruits ~ height*park, data=navsqu_all_sub, family=poisson(link="log"))

navsqu_all$fruits_pred <- ceiling(predict(m1, navsqu_all, type = "response"))
head(navsqu_all, n=8, width=Inf)
navsqu_all
ggplot(navsqu_all, aes(x=fruits, y=fruits_pred))+
  geom_point() + geom_smooth(method='lm')

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
navsqu_all <- navsqu_all %>%   # checked and don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))


## alone navsqu ----

navsqu_alone <- filter(alones_all, species=="NAVSQU") %>% 
  rowwise() %>%
  mutate(total.fr = rowSums(across(fr1:fr5), na.rm=TRUE),
         stems = ifelse(is.na(stems) & total.fr > 0, 1, stems),
         side_infl = max((total.fr - stems), 0),
         top_infl = min(stems, total.fr),
         fruits=topfr*top_infl + sidefr*side_infl)%>%
  ungroup()%>%
  as.data.frame()
head(navsqu_alone, n=3, width=Inf)

ggplot(navsqu_alone, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="navsqu alones")
ggsave("Figures/navsqu_ALONES_height_fruits_byPark.pdf", width=6, height=6)

ggplot(navsqu_alone, aes(x=height, y=fruits, color=park))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="navsqu alones")
ggsave("Figures/navsqu_ALONES_height_fruits.pdf", width=6, height=6)


navsqu_alone_sub <- navsqu_alone %>% filter(!(is.na(fruits)))  # fit fruit ~ height relationship with those plants with both
head(navsqu_alone_sub, n=4, width=Inf); dim(navsqu_all_sub)

# do these regressions by park
m2 <- glm(fruits ~ height*park , data=navsqu_alone_sub, family=poisson(link="log"))
summary(m2)

navsqu_alone$fruits_pred <- predict(m2, navsqu_alone, type = "response")   # type = "response"    needed to get fruits, not log(fruits)

ggplot(navsqu_alone, aes(x=fruits, y=fruits_pred))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
head(navsqu_alone); dim(navsqu_alone)
navsqu_alone <- navsqu_alone %>%   # checked and don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))%>% 
  dplyr::select(park:height, biomass:notes, fruits) %>% 
  as.data.frame()
navsqu_all <- navsqu_all %>% select(-c(stems,fruits_pred))

head(navsqu_alone, n=4, width=Inf); dim(navsqu_alone)
head(navsqu_all, n=4, width=Inf); dim(navsqu_all)

# EPIDEN ----
epiden_parks <- filter(focals, species=="EPIDEN") %>% rowwise() %>%
  mutate(fruits=exactfr,
         inflor=NA) %>%
  ungroup()
print(arrange(epiden_parks,desc(fruits)), n=8, width=Inf); dim(epiden_parks)

## Riverfront ----

head(filter(rf_focals, species=="EPIDEN"))
epiden_rf <- filter(rf_focals, species=="EPIDEN")
epiden_rf <- epiden_rf %>% rowwise() %>%
  mutate(fruits=exactfr,
         inflor=NA) %>%
  ungroup()
print(arrange(epiden_rf, desc(fruits)), n=8, width=Inf)

# so epiden_rf and epiden_parks have fruit estimates for plants with fruit counts

## get fruit~height relationships ----
# need to get estimates for the plants without fruit counts based on height  
# (doing this separately for 'diverse' plot plants versus alone plants)

# Put together diverse epiden (parks + rf)
epiden_all <- rbind(dplyr::select(epiden_rf, park:stems, biomass:notes, fruits), dplyr::select(epiden_parks, park:stems, biomass:notes, fruits)) %>%
  as.data.frame()
head(epiden_all, n=8, width=Inf)

ggplot(epiden_all, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="epiden")
ggsave("Figures/epiden_height_fruits_park.pdf", width=6, height=6)

ggplot(epiden_all, aes(x=height, y=fruits, color=plot))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="epiden")
ggsave("Figures/epiden_height_fruits.pdf", width=6, height=6)

ggplot(epiden_all, aes(x=height, y=biomass, color=plot))+
  geom_point() + geom_smooth(method='lm')+
  labs(title="epiden")
ggsave("Figures/epiden_height_biomass.pdf", width=6, height=6)


epiden_all_sub <- epiden_all %>% filter(!(is.na(fruits)))  # is.na(height) &
head(epiden_all_sub, n=4, width=Inf); dim(epiden_all_sub)

m3 <- glm(fruits ~ height*park, data=epiden_all_sub, family=poisson(link="log"))

epiden_all$fruits_pred <- ceiling(predict(m1, epiden_all, type = "response"))
head(epiden_all, n=8, width=Inf)
epiden_all

ggplot(epiden_all, aes(x=fruits_pred))+
  geom_histogram()
# one big outlier prediction for tall plant
epiden_all <- epiden_all %>%  mutate(fruits_pred = ifelse(fruits_pred>250, max(fruits)*2, fruits_pred))
ggplot(epiden_all, aes(x=fruits, y=fruits_pred))+
  geom_point() + geom_smooth(method='lm') + ylim(0,200)

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
epiden_all <- epiden_all %>%   # don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))
ggplot(epiden_all, aes(x=fruits))+
  geom_histogram()


## alone epiden ----

epiden_alone <- filter(alones_all, species=="EPIDEN") %>% 
  rowwise() %>%
  mutate(fruits=exactfr,
         inflor=NA)%>%
  ungroup()%>%
  as.data.frame()
head(epiden_alone, n=3, width=Inf)

ggplot(epiden_alone, aes(x=height, y=fruits, color=plot))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="epiden alones")
ggsave("Figures/epiden_ALONES_height_fruits_byPark.pdf", width=6, height=6)

ggplot(epiden_alone, aes(x=height, y=fruits, color=park))+
  # facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')+
  labs(title="epiden alones")
ggsave("Figures/epiden_ALONES_height_fruits.pdf", width=6, height=6)


epiden_alone_sub <- epiden_alone %>% filter(!(is.na(fruits)))  # fit fruit ~ height relationship with those plants with both
head(epiden_alone_sub, n=4, width=Inf); dim(epiden_all_sub)

# do these regressions by park
m2 <- glm(fruits ~ height*park , data=epiden_alone_sub, family=poisson(link="log"))
summary(m2)

epiden_alone$fruits_pred <- predict(m2, epiden_alone, type = "response")   # type = "response"    needed to get fruits, not log(fruits)

ggplot(epiden_alone, aes(x=fruits, y=fruits_pred))+
  facet_wrap(~park, scales = 'free')+
  geom_point() + geom_smooth(method='lm')

# keep the measured fruit counts and just use the predicted counts for plants that only have height?
head(epiden_alone); dim(epiden_alone)
epiden_alone <- epiden_alone %>%   # checked and don't need "rowwise" - works either way
  mutate(fruits = ifelse(is.na(fruits), fruits_pred, fruits))%>% 
  dplyr::select(park:height, biomass:notes, fruits) %>% 
  as.data.frame()
epiden_all <- epiden_all %>% select(-c(stems,fruits_pred))

head(epiden_alone, n=4, width=Inf); dim(epiden_alone)
head(epiden_all, n=4, width=Inf); dim(epiden_all)


# PLECON, GILCAP, CLAPUR, COLLIN, PLAFIG ----
# For PLECON, GILCAP, :
# Riverfront: FR1 is total inflorescences on the plant 
# Parks: total.fr is total inflorescences on the plant 
# Therefore, also need: fruits/infl --- get prior data

# For CLAPUR, COLLIN, PLAFIG:
# Riverfront: FR1 is total FRUITS on the plant 
# Parks: total.fr is total FRUITS on the plant 

# For EPIDEN:
# Riverfront: exactfr is total FRUITS on the plant 
# Parks: exactfr is total FRUITS on the plant 

# Seeds / fruit

sp2 <- c("PLECON", "CLAPUR", "COLLIN", "GILCAP", "PLAFIG")

head(focals,2); dim(focals)

# !! Placeholder estimate of 10 fruits / infl for PLEC & GILCAP ----
focals <- filter(focals, species %in% sp2)  %>% rowwise() %>% mutate(
  fruits = ifelse(species %in% sp2, total.fr, NA),
  inflor = NA,
  inflor = ifelse(species == "PLECON", total.fr, inflor),  # PLECON and GILCAP are just infl;  need fruits / infl
  inflor = ifelse(species == "GILCAP", total.fr, inflor),
  fruits = ifelse(species == "PLECON", inflor * 10, fruits),  # PLECON and GILCAP are just infl;  need fruits / infl
  fruits = ifelse(species == "GILCAP", inflor * 10, fruits)
) %>% ungroup() %>% as.data.frame()
head(focals,2); dim(focals)

rf_focals <- filter(rf_focals, species %in% sp2) %>% rowwise() %>% mutate(
  fruits = ifelse(species %in% sp2, fr1, NA),
  inflor = NA,
  inflor = ifelse(species == "PLECON", fr1, inflor),  # PLECON and GILCAP are just infl;  need fruits / infl
  inflor = ifelse(species == "GILCAP", fr1, inflor),
  fruits = ifelse(species == "PLECON", inflor * 10, fruits),  # PLECON and GILCAP are just infl;  need fruits / infl
  fruits = ifelse(species == "GILCAP", inflor * 10, fruits)) %>%
  ungroup() %>% as.data.frame()
head(rf_focals,2); dim(rf_focals)


# Alones
head(alones_all,2); dim(alones_all)
alones <- filter(alones_all, species %in% sp2) %>% rowwise() %>%
  mutate(fruits = ifelse(species %in% sp2, fr1, NA),
         inflor=NA,
         inflor = ifelse(species == "PLECON", fr1, inflor),  # PLECON and GILCAP are just infl;  need fruits / infl
         inflor = ifelse(species == "GILCAP", fr1, inflor),
         fruits = ifelse(species == "PLECON", inflor * 10, fruits),  # PLECON and GILCAP are just infl;  need fruits / infl
         fruits = ifelse(species == "GILCAP", inflor * 10, fruits)) %>% 
  ungroup() %>%   as.data.frame()
#
head(alones)

alone_summary <- alones %>% group_by(park, species) %>% 
  summarize(n_alones = length(which(fruits>0))) %>% ungroup()
pivot_wider(alone_summary, names_from=species, values_from=n_alones)
# No alones for COLLIN at SEM  (and we filtered them out because they were "GONE")

# Temp fix: simulate some 'alone' COLLIN plants at SEM with same mean fruits as WIR and BR
collin.temp <- filter(alones, species=="COLLIN", park %in% c("BR", "WIR")) 
collin.mean.frt.alone <- mean(collin.temp$fruits, na.rm=TRUE)
temp <- data.frame(
  park=rep("SEM",5), plot=rep("ALONE",5), 
  species=rep("COLLIN",5), rep=seq(1,5),date=rep(NA,5), plot.1=rep(NA,5), subplot=rep(NA,5), height=rep(NA,5), stems=rep(NA,5), fr1=rep(NA,5), fr2=rep(NA,5), fr3=rep(NA,5), fr4=rep(NA,5), fr5=rep(NA,5), exactfr=rep(NA,5), biomass=rep(NA,5), topfr=rep(NA,5), sidefr=rep(NA,5), notes=rep(NA,5), fruits=rpois(5,collin.mean.frt.alone), inflor=NA)

alones <- rbind(alones, temp)   
alone_summary <- alones %>% group_by(park, species) %>% 
  summarize(n_alones = length(which(fruits>0))) %>% ungroup()
pivot_wider(alone_summary, names_from=species, values_from=n_alones)

# ! why not GILCAP and PLECON in diverse plots at parks other than RF? ----
filter(focals, species=="PLECON", park=="SEM")
filter(focals, species=="GILCAP", park=="SEM")

# Merge what we have ----
head(collom_all, n=2, width=Inf); dim(collom_all)
head(navsqu_all, n=2, width=Inf); dim(navsqu_all)
head(focals, n=2)
head(rf_focals, n=2)
head(collom_alone, n=2, width=Inf); dim(collom_alone)
head(navsqu_alone, n=2, width=Inf); dim(navsqu_alone)
head(alones, n=2)
head(epiden_all, n=2, width=Inf); dim(epiden_all)
head(epiden_alone, n=2, width=Inf); dim(epiden_alone)

str(collom_all)
str(navsqu_all)
str(focals)
str(rf_focals)
str(collom_alone)
str(navsqu_alone)
str(alones)
str(epiden_all)
str(epiden_alone)

# columns to use
# park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits
dat_all <- rbind(
  dplyr::select(collom_all, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(navsqu_all, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(focals, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(rf_focals,  park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(collom_alone, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(navsqu_alone, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(alones, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(epiden_all, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits),
  dplyr::select(epiden_alone, park,plot,species,rep,date,plot.1,subplot,height,biomass,topfr,sidefr,notes,fruits)
)
head(dat_all); dim(dat_all)
str(dat_all)

# Assign seeds for each species ----
# !! Placeholder estimates of seeds: ----
# NAVSQU: 20 seeds / fruit 
# NAVSQU: 20 seeds / fruit 

no.fruits.data <- filter( dat_all, is.na(fruits))
head(no.fruits.data); dim(no.fruits.data)
## !! Work on getting values for these missing fruits   e.g. BT plants ----

head(dat_all); dim(dat_all)
dat_all <- filter(dat_all, !is.na(fruits))
head(dat_all); dim(dat_all)
## calculate seeds ----
dat_all <- dat_all %>% mutate(
  seeds = case_when(species=="COLLOM" ~ fruits * 3,
                    species=="NAVSQU" ~ fruits * 20, # placeholder
                    species=="PLECON" ~ fruits * 1,
                    species=="GILCAP" ~ fruits * 2,
                    species=="CLAPUR" ~ fruits * 20, # placeholder
                    species=="COLLIN" ~ fruits * 4,
                    species=="PLAFIG" ~ fruits * 4,
                    species=="EPIDEN" ~ fruits * 15) # placeholder
)
head(dat_all)

# Merge with Neighborhoods----
head(sem_neighborhoods,2); dim(sem_neighborhoods)
head(rf_neighborhoods,2); dim(rf_neighborhoods)
head(br_neighborhoods,2); dim(br_neighborhoods)
head(wir_neighborhoods,2); dim(wir_neighborhoods)

sem_neighborhoods <- sem_neighborhoods %>% mutate(unique=paste(park, plot, plot.1, subplot, sep='-'))
rf_neighborhoods <- rf_neighborhoods %>% mutate(unique=paste(park, plot, plot.1, subplot, sep='-'))
br_neighborhoods <- br_neighborhoods %>% mutate(unique=paste(park, plot, plot.1, subplot, sep='-'))
wir_neighborhoods <- wir_neighborhoods %>% mutate(unique=paste(park, plot, plot.1, subplot, sep='-'))

dat_all <- dat_all %>% mutate(unique=paste(park, plot, plot.1, subplot, sep='-'))

?dplyr::bind_rows
sem_dat <- left_join(filter(dat_all, park=="SEM"), sem_neighborhoods) 
rf_dat <- left_join(filter(dat_all, park=="RF"), rf_neighborhoods) 
br_dat <- left_join(filter(dat_all, park=="BR"), br_neighborhoods) 
wir_dat <- left_join(filter(dat_all, park=="WIR"), wir_neighborhoods) 

str(sem_dat)
str(rf_dat)
str(br_dat)
str(wir_dat)

filter(dat_all, unique=="RF-DIVERSE-1-C3")
filter(rf_dat, unique=="RF-DIVERSE-1-C3")
filter(rf_neighborhoods, unique=="RF-DIVERSE-1-C3")

# error: all the columns need to be the same class across the merge;  
# PLECON was chr in WIR  and integer in others; so removed "skip" from the data file


# how to join dataframes with overlapping but different columns
dat <- dplyr::bind_rows(sem_dat, rf_dat, br_dat, wir_dat) 
filter(dat, unique=="RF-DIVERSE-1-C3")
dim(dat)
table(dat$park)
head(filter(dat, plot=="ALONE"))

dat <- dat %>% mutate(N_neigh = rowSums(across(PLECON:ACMAME), na.rm=TRUE),
                      N_neigh8 = rowSums(across(PLECON:COLLOM), na.rm=TRUE)
                      )
# Plot raw curves ----

# Plot overall

ggplot(dat, aes(x=seeds, color=park)) +          
  geom_density() + xlim(0,1000)+ 
  facet_wrap(~ species , scales = "free") 
  
ggplot(dat, aes(x=N_neigh, y=log(seeds), color=park)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species , scales = "free") +
  # scale_color_manual(values = WA)+
  stat_smooth(method = "nls",
  formula = y ~ a/(1+b*x),
  method.args = list(start = list(a = 1000, b = .01)),
  se = FALSE) #+
ggsave("Figures/curves_overall_LOG.pdf", width=12, height=12)
# ggsave("Figures/curves_overall.pdf", width=12, height=12)

ggplot(dat, aes(x=N_neigh8, y=log(seeds), color=park)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species , scales = "free") +
  # scale_color_manual(values = WA)+
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) #+
ggsave("Figures/curves_overall_focalsOnly_LOG.pdf", width=12, height=12)
# ggsave("Figures/curves_overall_focalsOnly.pdf", width=12, height=12)

ggplot(dat, aes(x=N_neigh8, y=log(seeds), color=park)) +          
  geom_jitter(width=.1) + 
  facet_wrap(~ species , scales = "free") +
  xlim(0,50)+
  # scale_color_manual(values = WA)+
  stat_smooth(method = "nls",
              formula = y ~ a/(1+b*x),
              method.args = list(start = list(a = 1000, b = .01)),
              se = FALSE) #+
ggsave("Figures/curves_overall_focalsOnly_xlim50_LOG.pdf", width=12, height=12)
# ggsave("Figures/curves_overall_focalsOnly_xlim50.pdf", width=12, height=12)

# !! Apparant issues ----

# huge NAVSQU seeds value somewhere?
head(arrange(filter(dat, species=="NAVSQU"), desc(seeds)))


# Save data files ----

write.table(dat, "data_cleaned/dat.csv", sep=',', row.names=FALSE, eol = "\r")
write.table(sem_dat, "data_cleaned/sem_dat.csv", sep=',', row.names=FALSE, eol = "\r")
write.table(wir_dat, "data_cleaned/wir_dat.csv", sep=',', row.names=FALSE, eol = "\r")
write.table(br_dat, "data_cleaned/br_dat.csv", sep=',', row.names=FALSE, eol = "\r")
write.table(rf_dat, "data_cleaned/rf_dat.csv", sep=',', row.names=FALSE, eol = "\r")


