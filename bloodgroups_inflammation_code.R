#########################################
# SARAH VAN ALSTEN                      #
# NCI/DCEG; IIB                         #
# Date Created: July 21, 2020           #
# Purpose: Read in MIP Data and recode  #
# Packages used: tidyverse, haven,      #
# janitor, goft, corrplot, tableone, AER#
# Last Update: Sept 10, 2020            #
#########################################

#open packages needed (others we'll call a function or two explicitly with ::)
library(tidyverse)
library(janitor)


# Read in Data ------------------------------------------------------------
#note: in your working directory, have a folder called "data" which
#holds all of the MIP data files. Then these lines will work.

#read in the data: some are CSV's, some are sas7bdat
csvs <- paste0("data\\", list.files("data\\")[str_detect(list.files("data\\"), ".csv")])
sass <- paste0("data\\", list.files("data\\")[str_detect(list.files("data\\"), ".sas7bdat")])

csvs <- map(.x = csvs, .f = read.csv)
sass <- map(.x = sass, .f = haven::read_sas)

#add an identifier to state which study each is from
addID <- function(x, id){
  return(x %>% mutate(overall_study = id))
}

csvs <- map2(.x = csvs, .y = 1:8, .f = addID)
sass <- map2(.x = sass, .y = 9:10, .f = addID)

#extract individual studies from lists
lr <- csvs[[1]] #lung replication = 1
p2 <- csvs[[2]] #2nd pilot = 2
en <- csvs[[3]] #endometrial = 3
lu <- csvs[[4]] #first lung = 4
nh <- csvs[[5]] #nhl = 5
ov <- csvs[[6]] #ovarian = 6
p1 <- csvs[[7]] #1st pilot = 7
gi <- csvs[[8]] #upper gi = 8
co <- sass[[1]] #colorectal = 9
di <- sass[[2]] #diet qxn = 10

#sas files named things differently
names(co)[10] <- "X_6CKINE_OutsideLimits"


# Merge and Recode --------------------------------------------------------



#Bind rows for all of these together into 1 data set
#the mutates are to make sure that variable types are compatible
#warning lets us know that '.' got assigned NA as it should
dat <- bind_rows(lr %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric)) %>% clean_names(),
                 p2 %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 en %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 nh %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 ov %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 p1 %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 co %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 gi %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names(),
                 lu %>% mutate(across(contains("OutsideLimits"), .fns = as.numeric)) %>%
                   mutate(across(contains("pg_ml"), .fns = as.numeric))%>% clean_names())

#bind in the dietary/demographic info
dat <- dat %>% left_join(di, by = "plco_id")
dat <- dat %>% select(-overall_study.y) %>% rename(overall_study = overall_study.x)

#how many unique id ==> 4877
#length(unique(dat$plco_id))

#look at duplicate people
dupes <- get_dupes(dat, plco_id)
not.dupes <- dat %>% filter(!plco_id %in% dupes$plco_id)

#how many studies are each dupe in? ==> 220 in 2, 14 in 3
#table(dupes$dupe_count)

#which combos of studies?
#14 people with 3; 220 people with 2 = 234 total
#(219 of whom are white)
minDupe <- dupes %>%
  summarise(across(contains("pg_ml"), .fns = is.na)) %>%
  select_at(vars(contains("pg_ml")))

#add amount of missing data to the original duplicates
dupes$miss<- rowSums(minDupe[])

#select studies which have the LEAST missing data for each participant
dupes <- dupes %>%
  group_by(plco_id) %>%
  filter(miss == min(miss))

#that helps in some instances; others we'll pick the first study of inclusion
#for this, we need to sort by study date: lung-> NHL-> OV -> LR -> EN -> CO -> GI 
dupes <- dupes %>%
  mutate(study_order = case_when(overall_study == 4 ~ 1,
                                 overall_study == 5 ~ 2,
                                 overall_study == 6 ~ 3,
                                 overall_study == 1 ~ 4,
                                 overall_study == 3 ~ 5,
                                 overall_study == 9 ~ 6,
                                 overall_study == 8 ~ 7,
                                 TRUE ~ 8))
dupes <- dupes %>%
  arrange(study_order) %>%
  mutate(keep = lead(overall_study))

#non-duplicates are those where missing info was least
dupes.dupes <- get_dupes(dupes, plco_id)
dupes.nondupes <- dupes %>% filter(!plco_id %in% dupes.dupes$plco_id)

#for the ones where we look for first study:
dupes.dupes <- dupes.dupes %>%
  filter(!is.na(keep))

#get rid of all the extra columns I added so I can put back with 
#non-duplicates in original data
dupes.dupes <- dupes.dupes %>%
  select(-keep, -dupe_count, -miss, -study_order)
dupes.nondupes <- dupes.nondupes %>%
  select(-keep, -dupe_count, -miss, -study_order)


#now add back with non-duplicates into one dataframe
dat <- bind_rows(not.dupes, dupes.dupes, dupes.nondupes)

#some names of markers had different underscore patterns:
#make sure all of those are in ONE column
dat <- dat %>%
  mutate(sil6r_outside_limits = ifelse(!is.na(s_il_6r_outside_limits), 
                                       s_il_6r_outside_limits, sil_6r_outside_limits),
         sil4r_outside_limits = ifelse(!is.na(s_il_4r_outside_limits), 
                                       s_il_4r_outside_limits, sil_4r_outside_limits),
         silrii_outside_limits = ifelse(!is.na(s_ilrii_outside_limits), 
                                        s_ilrii_outside_limits, silrii_outside_limits),
         stnfrii_outside_limits = ifelse(!is.na(s_tnfrii_outside_limits), 
                                         s_tnfrii_outside_limits, stnfrii_outside_limits),
         stnfri_outside_limits = ifelse(!is.na(s_tnfri_outside_limits), 
                                        s_tnfri_outside_limits, stnfri_outside_limits),
         svegfr2_outside_limits = ifelse(!is.na(s_vegfr2_outside_limits), 
                                         s_vegfr2_outside_limits, svegfr2_outside_limits),
         svegfr3_outside_limits = ifelse(!is.na(s_vegfr3_outside_limits), 
                                         s_vegfr3_outside_limits, svegfr3_outside_limits),
         segfr_outside_limits = ifelse(!is.na(s_egfr_outside_limits), 
                                       s_egfr_outside_limits, segfr_outside_limits),
         tnfa_outside_limits = ifelse(!is.na(tn_fa_outside_limits), 
                                      tn_fa_outside_limits, tnfa_outside_limits),
         il_29_ifnl1_outside_limits = ifelse(!is.na(il_29_if_nl1_outside_limits), 
                                             il_29_if_nl1_outside_limits, il_29_ifnl1_outside_limits),
         sil6r_pg_ml = ifelse(!is.na(s_il_6r_pg_ml), 
                              s_il_6r_pg_ml, sil_6r_pg_ml),
         sil4r_pg_ml = ifelse(!is.na(s_il_4r_pg_ml), 
                              s_il_4r_pg_ml, sil_4r_pg_ml),
         silrii_pg_ml = ifelse(!is.na(s_ilrii_pg_ml), 
                               s_ilrii_pg_ml, silrii_pg_ml),
         stnfrii_pg_ml = ifelse(!is.na(s_tnfrii_pg_ml), 
                                s_tnfrii_pg_ml, stnfrii_pg_ml),
         stnfri_pg_ml = ifelse(!is.na(s_tnfri_pg_ml), 
                               s_tnfri_pg_ml, stnfri_pg_ml),
         svegfr2_pg_ml = ifelse(!is.na(s_vegfr2_pg_ml), 
                                s_vegfr2_pg_ml, svegfr2_pg_ml),
         svegfr3_pg_ml = ifelse(!is.na(s_vegfr3_pg_ml), 
                                s_vegfr3_pg_ml, svegfr3_pg_ml),
         segfr_pg_ml = ifelse(!is.na(s_egfr_pg_ml), 
                              s_egfr_pg_ml, segfr_pg_ml),
         tnfa_pg_ml = ifelse(!is.na(tn_fa_pg_ml), 
                             tn_fa_pg_ml, tnfa_pg_ml),
         x6ckine_pg_ml = ifelse(!is.na(x_6ckine_pg_ml), 
                                x_6ckine_pg_ml, x6ckine_pg_ml),
         x6ckine_outside_limits = x_6ckine_outside_limits,
         il_29_ifnl1_pg_ml = ifelse(!is.na(il_29_if_nl1_pg_ml), 
                                    il_29_if_nl1_pg_ml, il_29_ifnl1_pg_ml))

#now remove the old cols that I don't want w/ bad underscore pattern
dat <- dat %>% 
  select(-x_6ckine_outside_limits, -x_6ckine_pg_ml,
         -tn_fa_outside_limits, -tn_fa_pg_ml,
         -s_egfr_pg_ml, -s_egfr_outside_limits, 
         -s_vegfr2_outside_limits, -s_vegfr2_pg_ml,
         -s_vegfr3_pg_ml, -s_vegfr3_outside_limits,
         -s_tnfri_outside_limits, -s_tnfri_pg_ml,
         -s_tnfrii_pg_ml, -s_tnfrii_outside_limits,
         -s_ilrii_outside_limits, -s_ilrii_pg_ml,
         -s_il_4r_pg_ml, -s_il_4r_outside_limits,
         -s_il_6r_outside_limits, -s_il_6r_pg_ml,
         -il_29_if_nl1_pg_ml, -il_29_if_nl1_outside_limits)


#move columns so limit detect and markers are grouped together
dat <- dat %>%
  relocate(names(dat[,str_detect(names(dat), "outside_limits")]), .after = "sample_yr")

dat <- dat %>%
  relocate(names(dat[,str_detect(names(dat), "pg_ml")]), .after = "x6ckine_outside_limits" )

#finally, sort the markers/llod to be in same order
marks <- dat %>% select_at(vars(contains("pg_ml")))
llod <- dat %>% select_at(vars(contains("outside_limit")))
other <- dat %>% select(-contains("pg_ml")) %>% select(-contains("outside_limit"))

#put in alphabetical order
marks <- marks[ , sort(names(marks))]
llod <- llod[, sort(names(llod))]

#finally, put everything back together as dataframe
dat <- bind_cols(other, marks, llod)

#forgot to remove these old underscore patterns; do that now
dat <- dat %>%
  select(-sil_4r_outside_limits, -sil_4r_pg_ml,
         -sil_6r_outside_limits, -sil_6r_pg_ml)


#because we likely are going to want to only use white participants,
#subset based on race
dat <- dat %>% 
  filter(race7 == 1 & hispanic_f !=1)

# Get LLODs ---------------------------------------------------------------

f1 <- function(x){
  return(sum(x == -1, na.rm = T))
}
f2 <- function(x){
  return(sum(x == -1, na.rm = T)/sum(!is.na(x)))
}


belowLLOD <- dat %>%
  summarise(across(contains("outside_limit"), .fns = list(f1,f2), .names = "{col}.{fn}")) %>%
  pivot_longer(cols = 1:206) %>%
  mutate(name = str_remove_all(name, "\\.1")) %>%
  mutate(name = str_remove_all(name, "\\.2")) %>% 
  group_by(name) %>% 
  mutate(marker_id = row_number()) %>%
  pivot_wider(names_from = name, values_from = value) 


belowLLOD<- belowLLOD %>% t() %>% as.data.frame()
names(belowLLOD) <- c("count", "percent")
belowLLOD$marker <- rownames(belowLLOD)

belowLLOD <- belowLLOD[! belowLLOD$marker== "marker_id",]

belowLLOD[belowLLOD$percent >.9, "marker"]
#6 markers: IL3, LIF, SCD30, SILRI, SRAGE, SVEGFR1 > .90


# Get ULODs ----------------------------------------------------------------
f1 <- function(x){
  return(sum(x == 1, na.rm = T))
}
f2 <- function(x){
  return(sum(x == 1, na.rm = T)/sum(!is.na(x)))
}

aboveULOD <- dat %>%
  summarise(across(contains("outside_limit"), .fns = list(f1,f2), .names = "{col}.{fn}")) %>%
  pivot_longer(cols = 1:206) %>%
  mutate(name = str_remove_all(name, "\\.1")) %>%
  mutate(name = str_remove_all(name, "\\.2")) %>% 
  group_by(name) %>% 
  mutate(marker_id = row_number()) %>%
  pivot_wider(names_from = name, values_from = value) 

aboveULOD<- aboveULOD %>% t() %>% as.data.frame()
names(aboveULOD) <- c("count", "percent")
aboveULOD$marker <- rownames(aboveULOD)

aboveULOD <- aboveULOD[! aboveULOD$marker== "marker_id",]

aboveULOD[aboveULOD$percent >.9, "marker"] #no marker to exclude
aboveULOD[aboveULOD$percent >.7, "marker"] #one marker SCD40L

#next: which studies assessed which markers
studyMiss <- cbind.data.frame(
  dat %>% group_by(overall_study) %>%
    summarise(across(67:169, .fns = ~(sum(is.na(.))))),
  dat %>% group_by(overall_study) %>% count() %>% select(n))

#overall_study column got put on twice, remove
studyMiss <- studyMiss[, -105]

#what percent is missing within that study
studyMiss <- studyMiss %>%
  mutate(across(2:104, .fns = ~(./n)))


studyMiss <- studyMiss %>% t() %>% as.data.frame()
keep <- rownames(studyMiss)

table(dat$overall_study)

studyMiss <- studyMiss %>% 
  mutate(lr = ifelse(V1 < .9, "LR", ""),
         p2 = ifelse(V2 < .9, "P2", ""),
         lu = ifelse(V4 < .9, "LU", ""),
         nh= ifelse(V5 < .9, "NH", ""),
         ov = ifelse(V6 < .9, "OV", ""),
         p1  = ifelse(V7 < .9, "P1", ""),
         gi  = ifelse(V8 < .9, "GI", ""),
         co= ifelse(V9 < .9, "CO", ""),
         en  = ifelse(V3 < .9, "EN", ""))

studyMiss <- studyMiss %>%
  mutate(all = paste(lr, p2, en, lu, nh, ov, p1, gi, co, sep = " ")) %>%
  mutate(all = str_remove_all(all, "   "))


studyMiss$marker <- keep
studyMiss$marker <- str_remove_all(studyMiss$marker, "_pg_ml")
studyMiss$marker <- toupper(studyMiss$marker)

#Limit analyses to markers where <90% are below LLOD
dat <- dat %>%
  select(-contains(str_remove_all(belowLLOD[belowLLOD$percent >.9, "marker"], "outside_limits")))
# 
# for (i in 2:nrow(studyMiss)){
#   print(noquote(studyMiss$marker[i]))
# }
# for (i in 2:nrow(studyMiss)){
#   print(noquote(studyMiss$all[i]))
# }
# 
# for (i in 1:103){
#   print(noquote(sprintf(belowLLOD$percent[i]*100, fmt = "%.2f")))
# }
# 
# for (i in 1:103){
#   print(noquote(sprintf(aboveULOD$percent[i]*100, fmt = "%.2f")))
#   
# }

out <- cbind.data.frame(noquote(studyMiss$marker[2:104]),
                        noquote(studyMiss$all[2:104]),
                        noquote(sprintf(belowLLOD$percent*100, fmt = "%.2f")),
                        noquote(sprintf(aboveULOD$percent*100, fmt = "%.2f")))

names(out) <- c("Marker", "Studies", "Below_LLOD", "Above_LLOD")

out$Studies <- str_replace_all(out$Studies, "	 P2NH", "P2 NH") %>%
  str_replace_all(pattern = "P2NH", "P2 NH") %>%
  str_replace_all(pattern = " P2  LU", "P2 LU") %>%
  str_replace_all(pattern = "	ENCO", "EN CO") %>%
  str_replace_all(pattern = "	  EN  GI", "EN GI") %>%
  str_replace_all(pattern = "ENCO", "EN CO") %>%
  trimws()

#number of categories
out$Categories <- ifelse(as.numeric(out$Below_LLOD) <= 50, 4,
                         ifelse(as.numeric(out$Below_LLOD) <= 75, 3,
                                ifelse(as.numeric(out$Below_LLOD) <= 90, 2, NA)))
out$Categories <- ifelse(out$Marker == "SCD40L", 2, out$Categories)
out <- out %>% filter(Marker != "M_CSF")

dat <- dat %>% select(-m_csf_outside_limits, -m_csf_pg_ml)

#excluding those measured in only one study for the kruskal wallis tests
kw <- dat %>%
  select(-eotaxin_3_pg_ml, -ghrelin_pg_ml, -i_309_pg_ml,
         -il_20_pg_ml, -il_21_pg_ml, -il_28a_pg_ml,
         -il_9_pg_ml, -xcl1_lympho_pg_ml)

#variability by study: kruskal wallis tests
map(.x = kw[,68:154], .f = kruskal.test, g= kw$overall_study)
#only resistin, il23, and il13 have same distributions btwn studies

# Make Quartiles ----------------------------------------------------------
dat$overall_study <- as.factor(dat$overall_study)

#some markers had 0 below limit
zeroBelow <- str_replace_all(belowLLOD[belowLLOD$percent == 0, "marker"], "outside_limits", "pg_ml")[-30]
cat4 <- str_replace_all(belowLLOD[belowLLOD$percent <= .25, "marker"], "outside_limits", "pg_ml")[-30]
cat4b <- str_replace_all(belowLLOD[belowLLOD$percent > .25 & belowLLOD$percent <= .5, "marker"], "outside_limits", "pg_ml")
cat3 <- str_replace_all(belowLLOD[belowLLOD$percent > .5 & belowLLOD$percent <= .75, "marker"], "outside_limits", "pg_ml")
cat2 <- str_replace_all(belowLLOD[belowLLOD$percent > .75 & belowLLOD$percent <= .9, "marker"], "outside_limits", "pg_ml")

#for the 4 category (and for 4b), need to make sure that LLOD is consistent across
#studies; some individual studies might not have obs with the LLOD

#function to recode the columns based on categorization scheme
recodeLR <- function(x, cat){
  if (cat == "2"){
    return(ifelse(x == min(x, na.rm = T), 1, ifelse(!is.na(x), 2, NA)))
  } else if (cat == "3"){
    llod <- min(x, na.rm = T)
    x_no_llod <- x[x != llod]
    return(ifelse(x == llod, 1,
                  ifelse(x <= median(x_no_llod, na.rm = T), 2,
                         ifelse(!is.na(x), 3, NA))))
  } else if (cat == "4b"){
    llod <- min(x, na.rm = T)
    x_no_llod <- x[x != llod & ! is.na(x)]
    x_no_na <- x[!is.na(x)]
    return(ifelse(x == llod, 1,
                  ifelse(x <= quantile(x_no_llod, .333, na.rm = T), 2,
                         ifelse(x <= quantile(x_no_llod, .666, na.rm = T), 3, 
                                ifelse(!is.na(x), 4, NA)))))
    
  } else if (cat == "4"){
    llod <- min(x, na.rm = T)
    x_no_llod <- x[x != llod & ! is.na(x)]
    x_no_na <- x[!is.na(x)]
    return(ifelse(x <= quantile(x_no_llod, .25, na.rm = T), 1,
                  ifelse(x <= quantile(x_no_llod, .5, na.rm = T), 2,
                         ifelse(x <= quantile(x_no_llod, .75, na.rm = T), 3, 
                                ifelse(!is.na(x), 4, NA)))))
  }
}

#select markers that need to be categorized into 4/3/2 etc.
#and do so, renaming the category columns with "_cat" at end
test <- dat %>%
  group_by(overall_study) %>%
  #across the picked columns, apply function
  mutate(across(.cols = contains(cat4), .fns = recodeLR, cat = "4")) %>%
  #rename with _cat
  rename_at(vars(contains(cat4)), .funs = ~str_replace_all(.,"_pg_ml", "_cat")) %>%
  #only keep the new columns
  select_at(str_replace_all(cat4, "_pg_ml", "_cat")) %>%
  ungroup() %>% select(-overall_study)

test2 <- dat %>%
  group_by(overall_study) %>%
  mutate(across(.cols = contains(cat4b), .fns = recodeLR, cat = "4b")) %>%
  rename_at(vars(contains(cat4b)), .funs = ~str_replace_all(.,"_pg_ml", "_cat")) %>%
  select_at(str_replace_all(cat4b, "_pg_ml", "_cat"))%>%
  ungroup() %>% select(-overall_study)

test3 <- dat %>%
  group_by(overall_study) %>%
  mutate(across(.cols = contains(cat3), .fns = recodeLR, cat = "3")) %>%
  rename_at(vars(contains(cat3)), .funs = ~str_replace_all(.,"_pg_ml", "_cat")) %>%
  select_at(str_replace_all(cat3, "_pg_ml", "_cat"))%>%
  ungroup() %>% select(-overall_study)

test4 <- dat %>%
  group_by(overall_study) %>%
  mutate(across(.cols = contains(cat2), .fns = recodeLR, cat = "2")) %>%
  rename_at(vars(contains(cat2)), .funs = ~str_replace_all(.,"_pg_ml", "_cat")) %>%
  select_at(str_replace_all(cat2, "_pg_ml", "_cat"))%>%
  ungroup() %>% select(-overall_study)

#add the categorized columns to the DF
dat <- bind_cols(dat, test, test2, test3, test4)  

# Get Marker Conc by Quartile and Study -----------------------------------

markerCatNum <- dat %>%
  summarise(across(260:355, .fns = ~(length(levels(as.factor(.)))))) %>%
  pivot_longer(cols = 1:96)

names(markerCatNum) <- c("name", "value")
markerCatNum <- markerCatNum %>% arrange(name)

# #print it out
# for (i in 1:nrow(markerCatNum)){
#   print(noquote(markerCatNum$value[i]))
# }

fourCat <- markerCatNum %>%
  filter(value == 4) %>%
  dplyr::select(name) 
fourCat <- c(fourCat$name, str_replace(fourCat$name, "_cat", "_pg_ml"))

threeCat <- markerCatNum %>%
  filter(value == 3) %>%
  dplyr::select(name)
threeCat <-  c(threeCat$name, str_replace(threeCat$name, "_cat", "_pg_ml"))

twoCat <- markerCatNum %>%
  filter(value == 2) %>%
  dplyr::select(name)
twoCat <-  c(twoCat$name, str_replace(twoCat$name, "_cat", "_pg_ml"))

#set up a blank DF to hold results (32 cols across * number of markers down)
markerConc4s <- as.data.frame(cbind(c(fourCat[(length(fourCat)/2+1):length(fourCat)]),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2)),
                                    rep("", (length(fourCat)/2))))

names(markerConc4s) <- c("Marker", paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 1),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 2),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 3),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 4),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 5),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 6),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 7),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 8),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 9))

#make it a function
add.four <- function(i, study, cats){
  
  temp <- 
    dat %>%
    filter(overall_study == study)%>%
    select_at(fourCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+(length(fourCat)/2)-1), .funs = c(min,max), na.rm =T)
  
  #number skips 2 bc we aren't keepting pilot2
  #this tells which column to enter dat into
  addCol <- ifelse(study == 1, 0,
                   ifelse(study == 2, 4,
                          ifelse(study == 3, 8,
                                 ifelse(study == 4, 12,
                                        ifelse(study == 5, 16, 
                                               ifelse(study == 6, 20,
                                                      ifelse(study == 7, 24,
                                                             ifelse(study == 8, 28, 32))))))))
  
  markerConc4s[i, 2+addCol] <<- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc4s[i, 3+addCol] <<- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  markerConc4s[i, 4+addCol] <<- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  markerConc4s[i, 5+addCol] <<- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
}

markerConc3s <- as.data.frame(cbind(c(threeCat[(length(threeCat)/2+1):length(threeCat)]),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2)),
                                    rep("", (length(threeCat)/2))))

names(markerConc3s) <- c("Marker", paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 1),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 2),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 3),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 4),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 5),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 6),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 7),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 8),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 9))

#make it a function
add.three <- function(i, study, cats){
  
  temp <- 
    dat %>%
    filter(overall_study == study)%>%
    select_at(threeCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+(length(threeCat)/2)-1), .funs = c(min,max), na.rm =T)
  
  #number skips 2 bc we aren't keepting pilot2
  #this tells which column to enter dat into
  addCol <- ifelse(study == 1, 0,
                   ifelse(study == 2, 4,
                          ifelse(study == 3, 8,
                                 ifelse(study == 4, 12,
                                        ifelse(study == 5, 16, 
                                               ifelse(study == 6, 20,
                                                      ifelse(study == 7, 24,
                                                             ifelse(study == 8, 28, 32))))))))
  
  markerConc3s[i, 2+addCol] <<- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc3s[i, 3+addCol] <<- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  markerConc3s[i, 4+addCol] <<- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  #markerConc3s[i, 5+addCol] <<- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
}


markerConc2s <- as.data.frame(cbind(c(twoCat[(length(twoCat)/2+1):length(twoCat)]),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2)),
                                    rep("", (length(twoCat)/2))))

names(markerConc2s) <- c("Marker", paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 1),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 2),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 3),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 4),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 5),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 6),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 7),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 8),
                         paste0(c("Q1", "Q2", "Q3", "Q4"), sep = "_", 9))

#make it a function
add.two <- function(i, study, cats){
  
  temp <- 
    dat %>%
    filter(overall_study == study)%>%
    select_at(twoCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+(length(twoCat)/2)-1), .funs = c(min,max), na.rm =T)
  
  #number skips 2 bc we aren't keepting pilot2
  #this tells which column to enter dat into
  addCol <- ifelse(study == 1, 0,
                   ifelse(study == 2, 4,
                          ifelse(study == 3, 8,
                                 ifelse(study == 4, 12,
                                        ifelse(study == 5, 16, 
                                               ifelse(study == 6, 20,
                                                      ifelse(study == 7, 24,
                                                             ifelse(study == 8, 28, 32))))))))
  
  markerConc2s[i, 2+addCol] <<- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc2s[i, 3+addCol] <<- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  #markerConc3s[i, 4+addCol] <<- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  #markerConc3s[i, 5+addCol] <<- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
}


for (i in 1:9){ #number of studies
  for (j in 1:65){ #number of markers w/ that categorization (eg length fourCat/2)
    add.four(j, i, 4)
  }
}

for (i in 1:9){# number studies
  for (j in 1:18){ #number of markers w/ that categorization (eg length threeCat/2)
    add.three(j, i, 3)
  }
}

for (i in 1:9){ #number of studies
  for (j in 1:13){ #number of markers w/ that categorization (eg length twoCat/2)
    add.two(j, i, 2)
  }
}

#bind all of the marker data together
markerConcAll <- rbind(markerConc4s, markerConc3s, markerConc2s)

#add blank rows for the markers that we excluded
markerConcAll <- rbind.data.frame(markerConcAll,
                                  c("IL_3", rep("", 36)),
                                  c("LIF", rep("", 36)),
                                  c("SCD30", rep("", 36)),
                                  c("SIL1RI", rep("", 36)),
                                  #c("SIL2RA", rep("", 36)), #not needed in white only subsample
                                  c("SRAGE", rep("", 36)),
                                  c("SVEGFR1", rep("", 36)))

markerConcAll <- markerConcAll %>% arrange(Marker)
markerConcAll$Marker <- toupper(str_remove_all(markerConcAll$Marker, "_pg_ml"))

markerConcAll <- bind_cols(out, markerConcAll)
markerConcAll <- markerConcAll[,-6]
names(markerConcAll)[1] <- "Marker"
markerConcAll[markerConcAll == "NA - NA" | markerConcAll == "Inf - -Inf"] <- "--"


#write.csv(markerConcAll, "all_concentrations_new.csv")

# Recode Demographics -----------------------------------------------------
#recode race: 1 = white, 2 = black, 3 = hispanic, 4 = asian, 5 = Pacific Islander, 6 = American Indian
#recode education: 1 = <HS, 2 = HS, 3= Post HS but not College Grad, 4 = College +
dat <- dat %>%
  mutate(race = ifelse(is.na(hispanic_f), race7,
                       ifelse(hispanic_f == 1, 3, race7)),
         educ = ifelse(educat <= 2, 1, 
                       ifelse(educat == 3, 2, 
                              ifelse(educat < 6, 3, 
                                     ifelse(!is.na(educat), 4, NA))))) %>%
  mutate(race = as.factor(race),
         gender = as.factor(gender),
         educ = as.factor(educ),
         cig_stat = as.factor(cig_stat))

#frequencies of DZ's
dzFreq <- dat %>%
  select_at(vars(contains("_f"))) %>%
  select(-weight_f, - height_f, -match_fiscal_rndyear_group) %>%
  summarise_all(.funs = table, useNA = "ifany")

#percentages of DZ's
dzFreq <- rbind(dzFreq,
                dat %>%
                  select_at(vars(contains("_f"))) %>%
                  select(-weight_f, - height_f,  -match_fiscal_rndyear_group) %>%
                  summarise_all(.funs = table, useNA = "ifany") %>%
                  mutate_all(.funs = ~(. / nrow(dat))))

rownames(dzFreq) <- c("No", "Yes", "Missing", "NoPerc", "YesPerc", "MissingPerc")

#convert grams of drinks into number of standard drinks (14g = 1)
#then recode into categories (none/<weekly; moderate (1 for F, 2 for M); heavy (>1 F, >2 M))
dat <- dat %>%
  mutate(drinks = DT_ALC_ALC_DRINKS/14) %>%
  mutate(drinks = round(drinks, digits = 0)) %>%
  mutate(alc_stat = ifelse(drinks == 0, 0, 
                           ifelse(drinks < 3 & gender == "M", 1,
                                  ifelse(drinks < 2 & gender == "F", 1,
                                         ifelse(!is.na(drinks), 2, NA)))))
dat <- dat %>%
  mutate(age5 = case_when(age < 60 ~ 1,
                          age < 65 ~ 2,
                          age < 70 ~ 3,
                          age < 75 ~ 4
  ))
dat$age5 <- as.factor(dat$age5)


#if we were to exclude the pilot study, how consistent would distributions be?
#first make new KW dataset
kw_no_pilot <- kw %>% filter(overall_study != 7 & overall_study !=2)
map(.x = kw_no_pilot[,68:154], .f = possibly(kruskal.test, NA_real_), g= kw_no_pilot$overall_study)
#resistin, mcp3, il1a, il11 are only ones that don't differ


# Overall Concentration ---------------------------------------------------

dat <- dat %>% ungroup()
#also add categorizations for if we didn't do it specific to the individual studies
test <- dat %>%
  #across the picked columns, apply function
  mutate(across(.cols = contains(cat4), .fns = recodeLR, cat = "4")) %>%
  #rename with _cat
  rename_at(vars(contains(cat4)), .funs = ~str_replace_all(.,"_pg_ml", "_allC")) %>%
  #only keep the new columns
  select_at(str_replace_all(cat4, "_pg_ml", "_allC"))

test2 <- dat %>%
  #across the picked columns, apply function
  mutate(across(.cols = contains(cat4b), .fns = recodeLR, cat = "4b")) %>%
  #rename with _cat
  rename_at(vars(contains(cat4b)), .funs = ~str_replace_all(.,"_pg_ml", "_allC")) %>%
  #only keep the new columns
  select_at(str_replace_all(cat4b, "_pg_ml", "_allC"))

test3 <- dat %>%
  #across the picked columns, apply function
  mutate(across(.cols = contains(cat3), .fns = recodeLR, cat = "3")) %>%
  #rename with _cat
  rename_at(vars(contains(cat3)), .funs = ~str_replace_all(.,"_pg_ml", "_allC")) %>%
  #only keep the new columns
  select_at(str_replace_all(cat3, "_pg_ml", "_allC"))

test4 <- dat %>%
  #across the picked columns, apply function
  mutate(across(.cols = contains(cat2), .fns = recodeLR, cat = "2")) %>%
  #rename with _cat
  rename_at(vars(contains(cat2)), .funs = ~str_replace_all(.,"_pg_ml", "_allC")) %>%
  #only keep the new columns
  select_at(str_replace_all(cat2, "_pg_ml", "_allC"))


#add these new ones to the original data
dat <- bind_cols(dat, test, test2, test3, test4)




markerConc4 <- as.data.frame(cbind(c(fourCat[66:130]),
                                   rep("", 65),
                                   rep("", 65),
                                   rep("", 65),
                                   rep("", 65)))

names(markerConc4) <- c("Marker", "Q1", "Q2", "Q3", "Q4")

cats<- lapply(list(fourCat, threeCat, twoCat), str_replace_all, pattern = "_cat", replacement = "_allC") %>%
  unlist()

fourCat <- cats[1:130]
threeCat <- cats[131:166]
twoCat <- cats[167:192]

for (i in 1:65){
  temp <- 
    dat %>%
    select_at(fourCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+64), .funs = c(min,max), na.rm =T)
  
  markerConc4[i, 2] <- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc4[i, 3] <- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  markerConc4[i, 4] <- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  markerConc4[i, 5] <- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
  
}

#do the same for 3 category markers
markerConc3 <- as.data.frame(cbind(c(threeCat[18:36]),
                                   rep("", 18),
                                   rep("", 18),
                                   rep("", 18),
                                   rep("", 18)))

names(markerConc3) <- c("Marker", "Q1", "Q2", "Q3", "Q4")

for (i in 1:18){
  temp <- 
    dat %>%
    select_at(threeCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+17), .funs = c(min,max), na.rm =T)
  
  markerConc3[i, 2] <- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc3[i, 3] <- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  markerConc3[i, 4] <- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  #markerConc3[i, 5] <- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
  
}

#finally for the two category markers
markerConc2 <- as.data.frame(cbind(c(twoCat[14:26]),
                                   rep("", 13),
                                   rep("", 13),
                                   rep("", 13),
                                   rep("", 13)))

names(markerConc2) <- c("Marker", "Q1", "Q2", "Q3", "Q4")

for (i in 1:13){
  temp <- 
    dat %>%
    select_at(twoCat) %>%
    group_by_at(i) %>%
    summarise_at(.vars = (i+12), .funs = c(min,max), na.rm =T) #the summand should be -1 of the for loop
  
  markerConc2[i, 2] <- paste0(sprintf(temp[1,2], fmt = "%.1f"), " - ", sprintf(temp[1,3], fmt = "%.1f"))
  markerConc2[i, 3] <- paste0(sprintf(temp[2,2], fmt = "%.1f"), " - ", sprintf(temp[2,3], fmt = "%.1f"))
  #markerConc2[i, 4] <- paste0(sprintf(temp[3,2], fmt = "%.1f"), " - ", sprintf(temp[3,3], fmt = "%.1f"))
  #markerConc2[i, 5] <- paste0(sprintf(temp[4,2], fmt = "%.1f"), " - ", sprintf(temp[4,3], fmt = "%.1f"))
  
}

#bind all of the marker data together
markerConc <- rbind.data.frame(markerConc4, markerConc3, markerConc2)

#add blank rows for the markers that we excluded
markerConc <- rbind.data.frame(markerConc,
                               c("IL_3", rep("", 4)),
                               c("LIF", rep("", 4)),
                               c("SCD30", rep("", 4)),
                               c("SIL1RI", rep("", 4)),
                               #c("SIL2RA", rep("", 4)),
                               c("SRAGE", rep("", 4)),
                               c("SVEGFRI", rep("", 4)))
markerConc <- markerConc %>% arrange(Marker)
keep <- markerConc$Marker
markerConc <- markerConc[markerConc$Marker != "tnf_b_allC",]
markerConc <- markerConc[,-1]

#to make one single table that has overall concentrations, then
#study specific quartiles along w % above and below LLOD, put this together w/
#markerConcAll object

markerConc <- bind_cols(markerConcAll[,1:5],
                        markerConc,
                        markerConcAll[,6:ncol(markerConcAll)])


# Examine Correlations Among Markers --------------------------------------
library(corrplot) #make pretty plot


#log transform markers
logdat <- dat %>%
  mutate(across(contains("pg_ml"), .fns = log10)) %>%
  dplyr::select(contains("pg_ml"))
logdat <- logdat[, order(names(logdat))]
logdat$overall_study <- dat$overall_study
logdat$overall_study <- factor(logdat$overall_study)


#function to calculate partial correlation using residuals of linear models
getPCOR <- function(m1, m2){
  
  d <- logdat[, names(logdat) %in% c(m1, m2, "overall_study")]
  d <- na.omit(d)
  
  fmla1 <- paste(m1, "~ overall_study")
  fmla2 <- paste(m2, "~ overall_study")
  
  mm1 <- lm(fmla1, data = d)
  res1 <- mm1$residuals
  mm2 <- lm(fmla2, data = d)
  res2 <- mm2$residuals
  return(cor(res1,res2))
  
}

markCross <- cross_df(.l = list(arg1= names(logdat)[1:96], arg2= names(logdat)[1:96]))

partial.corr <- map2(.x = markCross$arg1,
                     .y = markCross$arg2,
                     possibly(getPCOR, NA_real_))

partial.corr <- matrix(unlist(partial.corr), nrow = 96)
colnames(partial.corr) <- str_remove(names(logdat)[1:96], "_pg_ml")
rownames(partial.corr) <- str_remove(names(logdat)[1:96], "_pg_ml")

partial.corr[upper.tri(partial.corr)] <-  " "


#correlation matrix lower triangle only
corr.mat <- cor(logdat, use = "pairwise.complete.obs")
corr.mat.copy <- corr.mat
corr.mat.copy[is.na(corr.mat.copy)] <- 0

corrplot(corr.mat.copy, method = "square", type = "upper",
         sig.level = .05,insig = "blank", order = "FPC",
         tl.col = "white", tl.cex =.25)

#make nice looking triangle which can go to supplement
corr.mat[upper.tri(corr.mat, diag = FALSE)] <- "" 
corr.mat <- as.data.frame(corr.mat)


# Test Log-Normality ------------------------------------------------------
map(dat %>% dplyr::select(contains("_pg_ml")), goft::lnorm_test)
#almost none are actually log normal via test, although histograms
#shown below show it may not be too bad









# Descriptives ------------------------------------------------------------
library(tableone)

dat_labels <- dat %>%
  mutate(overall_study = case_when(overall_study == 1 ~ "LR",
                                   overall_study == 2 ~ "P2",
                                   overall_study == 3 ~ "EN",
                                   overall_study == 4 ~ "LU",
                                   overall_study == 5 ~ "NH",
                                   overall_study == 6 ~ "OV",
                                   overall_study == 7 ~ "P1",
                                   overall_study == 8 ~ "GI",
                                   overall_study == 9 ~ "CO"))
dat_labels <- dat_labels %>%
  mutate(race = case_when(race == 1 ~ "White",
                          race == 2 ~ "Black",
                          race == 3 ~ "Hispanic",
                          race == 4 ~ "Asian",
                          race == 5 ~ "Pacific Islander",
                          race == 6 ~ "American Indian"),
         Smoking = case_when(cig_stat == 0  ~ "Never",
                             cig_stat == 1 ~ "Former",
                             cig_stat == 2 ~ "Current"),
         Age = case_when(age5 == 1 ~"55 - 59",
                         age5 ==2 ~ "60 - 64",
                         age5 == 3 ~ "65 - 69",
                         age5 == 4 ~ "70 - 74"),
         Education = case_when(educ == 1 ~ "< High School",
                               educ == 2 ~ "High School Degree",
                               educ == 3 ~ "Some College",
                               educ == 4 ~ "College Degree +"))

dat_labels <- dat_labels %>%
  mutate(Gender = gender)

#write.csv(print(CreateTableOne(vars = c("Gender", "Education", "Smoking", "Age"),
#                               data = dat_labels, strata = "overall_study")),
#          quote = FALSE, file = "unweightedtabletwo.csv")

CreateTableOne(vars = c("Gender", "Education", "Smoking", "Age", "is_case"),
               data = dat_labels, factorVars = "is_case", strata = "overall_study")

# # Graphs/Histograms -------------------------------------------------------
# markerNames <- names(dat)[68:163]
# makeDensity <- function(i, name){
#   df <- dat_labels[, c(18, i)]
#   names(df) <- c("study", "marker")
#   ggplot(df, aes(x = marker, fill = factor(study))) + geom_density(alpha = .1) +
#     xlab(name)
# }
# g <- map2(68:163, markerNames, makeDensity)
# g
# 
# #make density plots w/ log transformed markers
# makeDensityLog <- function(i, name){
#   df <- dat_labels[, c(18, i)]
#   names(df) <- c("study", "marker")
#   df$marker <- log(df$marker, base = 10)
#   ggplot(df, aes(x = marker, fill = factor(study))) + geom_density(alpha = .1)+
#     xlab(name)
# }
# g <- map2(68:163, markerNames, makeDensityLog)
# g
# 
# #possibly problematic ones:
# #vegf, tgfa, svegfr3, mip1b, mcp4, il7, gcsf, egf = bimodal
# 
# #write log transformed distributions to PDF File
# pdf("logtransformed_markers_by_study.pdf")
# g
# dev.off()
# 
# 
# #if we standardize first:
# makeDensityStand <- function(i, name){
#   df <- dat_labels[, c(18, i)]
#   names(df) <- c("study", "marker")
#   df$marker <- scale(df$marker, center = T, scale = T)
#   ggplot(df, aes(x = marker, fill = factor(study))) + geom_density(alpha = .1) +
#     xlab(name)
#   
# }
# g <- map2(68:163, markerNames, makeDensityStand)
# g
# 
# #if we log transform then standardize
# makeDensityLogStand <- function(i, name){
#   df <- dat_labels[, c(18, i)]
#   names(df) <- c("study", "marker")
#   df$marker <- log(df$marker, base = 10)
#   df$marker <- scale(df$marker, center = T, scale = T)
#   ggplot(df, aes(x = marker, fill = factor(study))) + geom_density(alpha = .1) +
#     xlab(name)
#   
# }
# g <- map2(68:163, markerNames, makeDensityLogStand)
# g
# 
# #write log transformed distributions to PDF File
# pdf("logtransformed_then_standardized_markers_by_study.pdf")
# g
# dev.off()
# 

# Mean/SD of Log-Transformed Markers by Study -----------------------------

#first, log-transform the markers
dat <- dat %>%
  mutate(across(contains("_pg_ml"), .fns = log10))

#prior to z-scoring: summarize the mean/sd by study
mean_sd <- dat %>%
  group_by(overall_study) %>%
  summarise(across(contains("_pg_ml"), .fns = ~(paste0(sprintf(mean(., na.rm =T), fmt = "%.1f"), " (",
                                                       sprintf(sd(., na.rm =T), fmt = "%.1f"), ")"))))

mean_sd[mean_sd == "NaN (NA)"] <- ""
mean_sd <- mean_sd %>% t() %>% as.data.frame()
names(mean_sd) <- c("LungRep", "Pilot2", "Endometrial", "Lung", "NHL", "Ovarian", "Pilot1",
                    "UpperGI", "Colorectal")
mean_sd <- mean_sd[-1,]
mean_sd$Marker <- rownames(mean_sd)
mean_sd$Marker <- toupper(str_remove_all(mean_sd$Marker, "_pg_ml"))
mean_sd <- mean_sd %>% relocate(Marker, .before = LungRep)

#add blank rows for the markers that we excluded
mean_sd <- rbind.data.frame(mean_sd,
                            c("IL_3", rep("", 9)),
                            c("LIF", rep("", 9)),
                            c("SCD30", rep("", 9)),
                            c("SIL1RI", rep("", 9)),
                            #c("SIL2RA", rep("", 9)),
                            c("SRAGE", rep("", 9)),
                            c("SVEGFRI", rep("", 9)))
mean_sd <- mean_sd %>% arrange(Marker)

markerConc <- bind_cols(markerConcAll[,1:5],
                        markerConc,
                        markerConcAll[,6:ncol(markerConcAll)])

mean_sd <- cbind.data.frame(mean_sd$Marker,
                            markerConc[1:102,2:4],
                            mean_sd[,2:10])



# Batch Variable for CRC --------------------------------------------------

crc <- haven::read_sas("data\\package-plco-646-3\\Colorectal Multiplex Immune Marker Panel Study\\mips_colo_custom_nov18_d090120.sas7bdat")
crc <- crc %>%
  select(plco_id, batch,C6P3_LotNum, AD5P1_LotNum, C16P2_LotNum,
         C17P1_LotNum, SR9P1_LotNum, CVD3P2_LotNum)

dat <- dat %>%
  left_join(crc, by = "plco_id")


# Read in ALL PLCO data ---------------------------------------------------

plco <- haven::read_sas("data\\package-plco-646-3\\Cohort dataset\\plco_646_cohort_mar20_d082620.sas7bdat")


#recode some of needed variables
#limit to NHW like in our study
plco <- plco %>% filter(race7 == 1)

plco <- plco %>%
  mutate(age5 = case_when(age < 60 ~ 1,
                          age < 65 ~ 2,
                          age < 70 ~ 3,
                          age < 75 ~ 4))

plco.m.case <- plco %>% filter(gender == "M" &
                                 (nhl_mips_is_case |colo_mips_is_case |
                                    ovar_mips_is_case |lung_mips_is_case |
                                    first_mips_is_case | upgi_mips_is_case | 
                                    endo_mips_is_case |lung2_mips_is_case))
plco.f.case <- plco %>% filter(gender == "F" &
                                 (nhl_mips_is_case |colo_mips_is_case |
                                    ovar_mips_is_case |lung_mips_is_case |
                                    first_mips_is_case | upgi_mips_is_case | 
                                    endo_mips_is_case |lung2_mips_is_case))

plco.m.cont <- plco %>% 
  filter(gender == "M" &
           (nhl_mips_in_study |colo_mips_in_study |
              ovar_mips_in_study |lung_mips_in_study |
              upgi_mips_in_study | 
              endo_mips_in_study |lung2_mips_in_study)) %>%
  filter(!(nhl_mips_is_case &colo_mips_is_case &
             ovar_mips_is_case &lung_mips_is_case &
             first_mips_is_case & upgi_mips_is_case & 
             endo_mips_is_case &lung2_mips_is_case & first_mips_is_case))

plco.f.cont <- plco %>% 
  filter(gender == "F" &
           (nhl_mips_in_study |colo_mips_in_study |
              ovar_mips_in_study |lung_mips_in_study |
              upgi_mips_in_study | 
              endo_mips_in_study |lung2_mips_in_study)) %>%
  filter(!(nhl_mips_is_case &colo_mips_is_case &
             ovar_mips_is_case &lung_mips_is_case &
             first_mips_is_case & upgi_mips_is_case & 
             endo_mips_is_case &lung2_mips_is_case & first_mips_is_case))


#were they a control?
plco <- plco %>%
  mutate(colo_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id) & colo_mips_in_study, 1, 0),
         nhl_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& nhl_mips_in_study, 1, 0),
         upgi_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& upgi_mips_in_study, 1, 0),
         endo_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& endo_mips_in_study, 1, 0),
         lung_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& lung_mips_in_study, 1, 0),
         lung2_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& lung2_mips_in_study, 1, 0),
         ovar_mips_is_control = ifelse(plco_id %in% c(plco.f.cont$plco_id, plco.m.cont$plco_id)& ovar_mips_in_study, 1, 0))



# plco <- plco %>%
#   mutate(is_case = case_when(ovar_mips_is_case ==1 | lung2_mips_is_case ==1 | endo_mips_is_case == 1 |
#                                upgi_mips_is_case == 1 | lung_mips_is_case == 1 | nhl_mips_is_case == 1 |
#                                colo_mips_is_case == 1 ~ 1,
#                              T ~ 0))
# plco <- plco %>%
#   mutate(is_control = case_when(ovar_mips_is_control ==1 | lung2_mips_is_control ==1 | endo_mips_is_control == 1 |
#                                upgi_mips_is_control == 1 | lung_mips_is_control == 1 | nhl_mips_is_control == 1 |
#                                colo_mips_is_control == 1 ~ 1,
#                              T ~ 0))


had.draw <- plco %>% filter(had_blood_draw==1)

table(cont =had.draw$is_control, had.draw$is_case)
table(had.draw$is_case, had.draw$gender)

#FUNCTION TO GET IPW 
getIPW <- function(mod, variable){
  
  predictions <- predict(mod, type = "response")
  
  #define probs
  p_actual <- (mod$data[, names(mod$data) == variable] * predictions)+
    ((1-mod$data[, names(mod$data) == variable] )*(1-predictions))
  
  #define IPW
  ipw_case <- 1/p_actual
  
  ipw_case_stabilized <- ifelse(str_detect(variable, "case"), 
                                (mean(predictions))/predictions,
                                (mean(1-predictions))/(1-predictions))
  
  # TRIMMING the stabilized IPTW weights
  # dat$stable.trim.1.iptw <- ifelse(dat$stable.1.iptw < .10, 0.10, dat$stable.1.iptw)
  # dat$stable.trim.1.iptw <- ifelse(dat$stable.1.iptw > 10, 10, dat$stable.trim.1.iptw)
  # 
  # dat$stable.trim.2.iptw <- ifelse(dat$stable.2.iptw < .10, 0.10, dat$stable.2.iptw)
  # dat$stable.trim.2.iptw <- ifelse(dat$stable.2.iptw > 10, 10, dat$stable.trim.2.iptw)
  # 
  # dat$stable.trim.3.iptw <- ifelse(dat$stable.3.iptw < .10, 0.10, dat$stable.3.iptw)
  # dat$stable.trim.3.iptw <- ifelse(dat$stable.3.iptw > 10, 10, dat$stable.trim.3.iptw)
  
  
  ipw_case <- cbind.data.frame(ipw_case, ipw_case_stabilized, mod$data$plco_id)
  names(ipw_case) <- c("ipw", "ipw_stabilized", "plco_id")
  
  
  return(ipw_case)
}


#BEING IN THE CRC STUDY
#########################################################################
crc_mod_m_case <- glm(colo_mips_is_case ~ factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "M" & !plco$crohn_f & plco$has_adequate_fsg &
                                    !plco$colitis_f & !plco$polypos_f& !plco$gardner_f & plco$had_blood_draw,], 
                      family = "binomial")
crc_mod_f_case <- glm(colo_mips_is_case ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "F" & !plco$crohn_f & plco$has_adequate_fsg &
                                    !plco$colitis_f & !plco$polypos_f& !plco$gardner_f& plco$had_blood_draw,], 
                      family = "binomial")
crc_mod_m_cont <- glm(colo_mips_is_control ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "M" & !plco$crohn_f & plco$has_adequate_fsg &
                                    !plco$colitis_f & !plco$polypos_f& !plco$gardner_f& plco$had_blood_draw,], 
                      family = "binomial")
crc_mod_f_cont <- glm(colo_mips_is_control ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "F" & !plco$crohn_f & plco$has_adequate_fsg &
                                    !plco$colitis_f & !plco$polypos_f & !plco$gardner_f& plco$had_blood_draw,], 
                      family = "binomial")

ipw_colo_f_case <- getIPW(crc_mod_f_case, "colo_mips_is_case")%>% drop_na(ipw)
ipw_colo_f_cont <- getIPW(crc_mod_f_cont, "colo_mips_is_control")%>% drop_na(ipw)
ipw_colo_m_case <- getIPW(crc_mod_m_case, "colo_mips_is_case")%>% drop_na(ipw)
ipw_colo_m_cont <- getIPW(crc_mod_m_cont, "colo_mips_is_control")%>% drop_na(ipw)



#BEING IN OVARIAN STUDY
################################################################################
ova_mod_f_case <- glm(ovar_mips_is_case ~ factor(cig_stat) + factor(age5) +
                        pack_years, data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                      family = "binomial")
ova_mod_f_cont <- glm(ovar_mips_is_control ~ factor(cig_stat) + factor(age5) +
                        pack_years, data = plco[plco$gender == "F" & plco$had_blood_draw & (plco$ovariesr_f == 0),], 
                      family = "binomial")

ipw_ovar_f_case <- getIPW(ova_mod_f_case, "ovar_mips_is_case")%>% drop_na(ipw)
ipw_ovar_f_cont <- getIPW(ova_mod_f_cont, "ovar_mips_is_control")%>% drop_na(ipw)

#################################################################################
#BEING IN ENDOMETRIAL STUDY
end_mod_f_case <- glm(endo_mips_is_case ~ factor(cig_stat) + factor(age5) +
                        pack_years, data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                      family = "binomial")
end_mod_f_cont <- glm(endo_mips_is_control ~ factor(cig_stat) + factor(age5) +
                        pack_years, data = plco[plco$gender == "F"& plco$had_blood_draw & !plco$hyster_f,], 
                      family = "binomial")

ipw_endo_f_case <- getIPW(end_mod_f_case, "endo_mips_is_case")%>% drop_na(ipw)
ipw_endo_f_cont <- getIPW(end_mod_f_cont, "endo_mips_is_control")%>% drop_na(ipw)
##############################################################################
#BEING IN NHL STUDY
nhl_mod_m_case <- glm(nhl_mips_is_case ~ factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                      family = "binomial")
nhl_mod_f_case <- glm(nhl_mips_is_case ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                      family = "binomial")
nhl_mod_m_cont <- glm(nhl_mips_is_control ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                      family = "binomial")
nhl_mod_f_cont <- glm(nhl_mips_is_control ~  factor(smkstp) + factor(age5),
                      data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                      family = "binomial")

ipw_nhl_f_case <- getIPW(nhl_mod_f_case, "nhl_mips_is_case")%>% drop_na(ipw)
ipw_nhl_f_cont <- getIPW(nhl_mod_f_cont, "nhl_mips_is_control")%>% drop_na(ipw)
ipw_nhl_m_case <- getIPW(nhl_mod_m_case, "nhl_mips_is_case")%>% drop_na(ipw)
ipw_nhl_m_cont <- getIPW(nhl_mod_m_cont, "nhl_mips_is_control")%>% drop_na(ipw)
##############################################################################
#BEING IN LUNG STUDY
lung_mod_m_case <- glm(lung_mips_is_case ~ factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                       family = "binomial")
lung_mod_f_case <- glm(lung_mips_is_case ~  factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                       family = "binomial")
lung_mod_m_cont <- glm(lung_mips_is_control ~  factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                       family = "binomial")
lung_mod_f_cont <- glm(lung_mips_is_control ~  factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                       family = "binomial")

ipw_lung_f_case <- getIPW(lung_mod_f_case, "lung_mips_is_case")%>% drop_na(ipw)
ipw_lung_f_cont <- getIPW(lung_mod_f_cont, "lung_mips_is_control")%>% drop_na(ipw)
ipw_lung_m_case <- getIPW(lung_mod_m_case, "lung_mips_is_case")%>% drop_na(ipw)
ipw_lung_m_cont <- getIPW(lung_mod_m_cont, "lung_mips_is_control")%>% drop_na(ipw)

##############################################################################
#BEING IN lung2 STUDY
lung2_mod_m_case <- glm(lung2_mips_is_case ~ factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                        family = "binomial")
lung2_mod_f_case <- glm(lung2_mips_is_case ~  factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                        family = "binomial")
lung2_mod_m_cont <- glm(lung2_mips_is_control ~  factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                        family = "binomial")
lung2_mod_f_cont <- glm(lung2_mips_is_control ~  factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                        family = "binomial")

ipw_lung2_f_case <- getIPW(lung2_mod_f_case, "lung2_mips_is_case")%>% drop_na(ipw)
ipw_lung2_f_cont <- getIPW(lung2_mod_f_cont, "lung2_mips_is_control")%>% drop_na(ipw)
ipw_lung2_m_case <- getIPW(lung2_mod_m_case, "lung2_mips_is_case")%>% drop_na(ipw)
ipw_lung2_m_cont <- getIPW(lung2_mod_m_cont, "lung2_mips_is_control")%>% drop_na(ipw)

##############################################################################
#BEING IN UPGI STUDY
upgi_mod_m_case <- glm(upgi_mips_is_case ~ factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                       family = "binomial")
upgi_mod_f_case <- glm(upgi_mips_is_case ~  factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                       family = "binomial")
upgi_mod_m_cont <- glm(upgi_mips_is_control ~  factor(smkstp) + factor(age5),
                       data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                       family = "binomial")
upgi_mod_f_cont <- glm(upgi_mips_is_control ~ factor(age5),
                       data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                       family = "binomial")

ipw_upgi_f_case <- getIPW(upgi_mod_f_case, "upgi_mips_is_case")%>% drop_na(ipw)
ipw_upgi_f_cont <- getIPW(upgi_mod_f_cont, "upgi_mips_is_control")%>% drop_na(ipw)
ipw_upgi_m_case <- getIPW(upgi_mod_m_case, "upgi_mips_is_case")%>% drop_na(ipw)
ipw_upgi_m_cont <- getIPW(upgi_mod_m_cont, "upgi_mips_is_control")%>% drop_na(ipw)
##########################################################################
#BEING IN PILOT
first_mod_m_case <- glm(first_mips_is_case ~ factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "M"& plco$had_blood_draw,], 
                        family = "binomial")
first_mod_f_case <- glm(first_mips_is_case ~  factor(smkstp) + factor(age5),
                        data = plco[plco$gender == "F"& plco$had_blood_draw,], 
                        family = "binomial")

ipw_first_f_case <- getIPW(first_mod_f_case, "first_mips_is_case") %>% drop_na(ipw)
ipw_first_m_case <- getIPW(first_mod_m_case, "first_mips_is_case")%>% drop_na(ipw)


#join everything together: keeping ONLY the plco_ids that are already in data
dat <- dat %>%
  left_join(ipw_colo_f_case, by = "plco_id") %>%
  rename(colo_f_case_weight = ipw) %>%
  left_join(ipw_colo_f_cont, by = "plco_id") %>%
  rename(colo_f_cont_weight = ipw) %>%
  left_join(ipw_colo_m_cont, by = "plco_id") %>%
  rename(colo_m_cont_weight = ipw) %>%
  left_join(ipw_colo_m_case, by = "plco_id") %>%
  rename(colo_m_case_weight = ipw)

dat <- dat %>%
  left_join(ipw_lung_f_case, by = "plco_id") %>%
  rename(lung_f_case_weight = ipw) %>%
  left_join(ipw_lung_f_cont, by = "plco_id") %>%
  rename(lung_f_cont_weight = ipw) %>%
  left_join(ipw_lung_m_cont, by = "plco_id") %>%
  rename(lung_m_cont_weight = ipw) %>%
  left_join(ipw_lung_m_case, by = "plco_id") %>%
  rename(lung_m_case_weight = ipw)


dat <- dat %>%
  left_join(ipw_lung2_f_case, by = "plco_id") %>%
  rename(lung2_f_case_weight = ipw) %>%
  left_join(ipw_lung2_f_cont, by = "plco_id") %>%
  rename(lung2_f_cont_weight = ipw) %>%
  left_join(ipw_lung2_m_cont, by = "plco_id") %>%
  rename(lung2_m_cont_weight = ipw) %>%
  left_join(ipw_lung2_m_case, by = "plco_id") %>%
  rename(lung2_m_case_weight = ipw)

dat <- dat %>%
  left_join(ipw_nhl_f_case, by = "plco_id") %>%
  rename(nhl_f_case_weight = ipw) %>%
  left_join(ipw_nhl_f_cont, by = "plco_id") %>%
  rename(nhl_f_cont_weight = ipw) %>%
  left_join(ipw_nhl_m_cont, by = "plco_id") %>%
  rename(nhl_m_cont_weight = ipw) %>%
  left_join(ipw_nhl_m_case, by = "plco_id") %>%
  rename(nhl_m_case_weight = ipw)

dat <- dat %>%
  left_join(ipw_upgi_f_case, by = "plco_id") %>%
  rename(upgi_f_case_weight = ipw) %>%
  left_join(ipw_upgi_f_cont, by = "plco_id") %>%
  rename(upgi_f_cont_weight = ipw) %>%
  left_join(ipw_upgi_m_cont, by = "plco_id") %>%
  rename(upgi_m_cont_weight = ipw) %>%
  left_join(ipw_upgi_m_case, by = "plco_id") %>%
  rename(upgi_m_case_weight = ipw)

dat <- dat %>%
  left_join(ipw_endo_f_case, by = "plco_id") %>%
  rename(endo_f_case_weight = ipw) %>%
  left_join(ipw_endo_f_cont, by = "plco_id") %>%
  rename(endo_f_cont_weight = ipw)

dat <- dat %>%
  left_join(ipw_ovar_f_case, by = "plco_id") %>%
  rename(ovar_f_case_weight = ipw) %>%
  left_join(ipw_ovar_f_cont, by = "plco_id") %>%
  rename(ovar_f_cont_weight = ipw)

dat <- dat %>%
  left_join(ipw_first_f_case, by = "plco_id") %>%
  rename(first_f_case_weight = ipw) %>%
  left_join(ipw_first_m_case, by = "plco_id") %>%
  rename(first_m_case_weight = ipw)


#define single study-specific weights
dat <- dat %>%
  mutate(ipw_1 = case_when(gender == "M" & is_case ~lung2_m_case_weight,
                           gender == "M" & !is_case ~lung2_m_cont_weight,
                           gender == "F" & is_case ~lung2_f_case_weight,
                           gender == "F" & !is_case ~lung2_f_cont_weight),
         ipw_3 = case_when(
           gender == "F" & is_case ~endo_f_case_weight,
           gender == "F" & !is_case ~endo_f_cont_weight),
         ipw_4 = case_when(gender == "M" & is_case ~lung_m_case_weight,
                           gender == "M" & !is_case ~lung_m_cont_weight,
                           gender == "F" & is_case ~lung_f_case_weight,
                           gender == "F" & !is_case ~lung_f_cont_weight),
         ipw_5 = case_when(gender == "M" & is_case ~nhl_m_case_weight,
                           gender == "M" & !is_case ~nhl_m_cont_weight,
                           gender == "F" & is_case ~nhl_f_case_weight,
                           gender == "F" & !is_case ~nhl_f_cont_weight),
         ipw_6 = case_when(
           gender == "F" & is_case ~ovar_f_case_weight,
           gender == "F" & !is_case ~ovar_f_cont_weight),
         ipw_7 = case_when(gender == "M"  ~first_m_case_weight,
                           gender == "F"  ~first_f_case_weight),
         ipw_8 = case_when(gender == "M" & is_case ~upgi_m_case_weight,
                           gender == "M" & !is_case ~upgi_m_cont_weight,
                           gender == "F" & is_case ~upgi_f_case_weight,
                           gender == "F" & !is_case ~upgi_f_cont_weight),
         ipw_9 = case_when(gender == "M" & is_case ~colo_m_case_weight,
                           gender == "M" & !is_case ~colo_m_cont_weight,
                           gender == "F" & is_case ~colo_f_case_weight,
                           gender == "F" & !is_case ~colo_f_cont_weight))

#write.csv(dat, "data_with_ipw0910.csv")

#trim weights
dat <- dat %>%
  mutate(across(starts_with("ipw_"), .fns = ~(ifelse(. > 1000, 1000, .))))


#remove excess spaces from strings
markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "LR EN LU NH  GI CO",
                                         "LR EN LU NH GI CO")

markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "LUGI",
                                         "LU GI")

markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "LR EN LU  OV GI CO",
                                         "LR EN LU OV GI CO")

markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "EN  CO",
                                         "EN CO")
markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "LUP1",
                                         "LU P1")
markerConcAll$Studies <- str_replace_all(markerConcAll$Studies, "  ",
                                         " ")


unique(markerConcAll$Studies) #18 different types

#gather the subsets of markers that each collection corresponds to
m1 <- markerConcAll[markerConcAll$Studies == unique(markerConcAll$Studies)[1], "Marker"]
m2 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[2], "Marker"]
m3 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[3], "Marker"]
m4 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[4], "Marker"]
m5 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[5], "Marker"]
m6 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[6], "Marker"]
m7 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[7], "Marker"]
m8 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[8], "Marker"]
m9 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[9], "Marker"]
m10 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[10], "Marker"]
m11 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[11], "Marker"]
m12 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[12], "Marker"]
m13 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[13], "Marker"]
m14 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[14], "Marker"]
m15 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[15], "Marker"]
m16 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[16], "Marker"]
m17 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[17], "Marker"]
m18 <- markerConcAll[markerConcAll$Studies ==  unique(markerConcAll$Studies)[18], "Marker"]

m1 <- tolower(paste0(m1, "_pg_ml"))
m2 <- tolower(paste0(m2, "_pg_ml"))
m3 <- tolower(paste0(m3, "_pg_ml"))
m4 <- tolower(paste0(m4, "_pg_ml"))
m5 <- tolower(paste0(m5, "_pg_ml"))
m6 <- tolower(paste0(m6, "_pg_ml"))
m7 <- tolower(paste0(m7, "_pg_ml"))
m8 <- tolower(paste0(m8, "_pg_ml"))
m9 <- tolower(paste0(m9, "_pg_ml"))
m10 <- tolower(paste0(m10, "_pg_ml"))
m11 <- tolower(paste0(m11, "_pg_ml"))
m12 <- tolower(paste0(m12, "_pg_ml"))
m13 <- tolower(paste0(m13, "_pg_ml"))
m14 <- tolower(paste0(m14, "_pg_ml"))
m15 <- tolower(paste0(m15, "_pg_ml"))
m16 <- tolower(paste0(m16, "_pg_ml"))
m17 <- tolower(paste0(m17, "_pg_ml"))
m18 <- tolower(paste0(m18, "_pg_ml"))


#Marker specific IPW are generated as average of studies that the marker
#is included in (sep for males and females). Generate the 18 different IPWs for
#marker groupings


#m1 consists of EN, GI, CO. Avg of 3 for Female, avg of gi co for males
dat <- dat %>%
  mutate(m1_ipw = ifelse(gender == "F", (ipw_3+ipw_8+ipw_9)/3, (ipw_8+ipw_9)/2))

#m2 consists of P2, NH, OV, P1
dat <- dat %>%
  mutate(m2_ipw = ifelse(gender == "F", (ipw_6+ipw_5 + ipw_7)/3, (ipw_5  + ipw_7)/2))

#m3 consists of LR P2 EN LU NH P1 GI CO
dat <- dat %>%
  mutate(m3_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_7+ipw_8+ipw_9)/7,
                         (ipw_1+ipw_4+ipw_5+ipw_7+ipw_8+ipw_9)/6))

#m4 consists of LR P2 EN LU OV P1 GI CO
dat <- dat %>%
  mutate(m4_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_6+ipw_7+ipw_8+ipw_9)/7,
                         (ipw_1+ipw_4+ipw_7+ipw_8+ipw_9)/5))

#m5 LR P2 EN LU NH OV P1 GI CO
dat <- dat %>%
  mutate(m5_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_6+ipw_7+ipw_8+ipw_9)/8,
                         (ipw_1+ipw_4+ipw_5+ipw_7+ipw_8+ipw_9)/6))

#m6 P1
dat <- dat %>%
  mutate(m6_ipw = ipw_7)

#m7 P2 LU NH OV P1
dat <- dat %>%
  mutate(m7_ipw = ifelse(gender == "F", (ipw_5+ipw_6+ipw_7)/3,
                         (ipw_5+ipw_7)/2))

#m8 P2 LU NH OV P1 GI
dat <- dat %>%
  mutate(m8_ipw = ifelse(gender == "F", (ipw_5+ipw_6+ipw_7+ipw_8)/4,
                         (ipw_5+ipw_7+ipw_8)/3))

#m9 LU NH OV P1 GI
dat <- dat %>%
  mutate(m9_ipw = ifelse(gender == "F", (ipw_5+ipw_6+ipw_7+ipw_8)/4,
                         (ipw_5+ipw_7+ipw_8)/3))

#m10 P2 LU NH P1
dat <- dat %>%
  mutate(m10_ipw = (ipw_4+ipw_5+ipw_7)/3)


#m11 P1 GI
dat <- dat %>%
  mutate(m11_ipw = (ipw_7+ipw_9)/2)

#m12 LU NH OV P1
dat <- dat %>%
  mutate(m12_ipw = ifelse(gender == "F", (ipw_4+ipw_5+ipw_6+ipw_7)/4,
                          (ipw_4+ipw_5+ipw_7)/3))

#m13 LR EN LU NH OV P1 GI CO
dat <- dat %>%
  mutate(m13_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_6+ipw_7+ipw_8+ipw_9)/8,
                          (ipw_1+ipw_4+ipw_5+ipw_7+ipw_8+ipw_9)/6))

#m14 GI
dat <- dat %>%
  mutate(m14_ipw = ipw_8)

#m15 LR EN LU NH OV GI CO
dat <- dat %>%
  mutate(m15_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_6+ipw_8+ipw_9)/7,
                          (ipw_1+ipw_4+ipw_5+ipw_8+ipw_9)/5))

#m16 LU NH OV GI
dat <- dat %>%
  mutate(m16_ipw = ifelse(gender == "F", (ipw_4+ipw_5+ipw_6+ipw_8)/4,
                          (ipw_4+ipw_5+ipw_8)/3))

#m17 LR P2 EN LU P1 GI CO
dat <- dat %>%
  mutate(m17_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_7+ipw_8+ipw_9)/6,
                          (ipw_1+ipw_4+ipw_7+ipw_8+ipw_9)/5))

#m18 EN CO
dat <- dat %>%
  mutate(m18_ipw = ifelse(gender == "F", (ipw_3+ipw_9)/2, ipw_9))

#write.csv(dat, "analysis_data_with_constructed_weights.csv")

#dat <- read.csv("final_analysis_data_0912.csv")

#sum of weights
sum(dat$ipw_1, na.rm =T)
sum(dat$ipw_2, na.rm =T)
sum(dat$ipw_3, na.rm =T)
sum(dat$ipw_4, na.rm =T)
sum(dat$ipw_5, na.rm =T)
sum(dat$ipw_6, na.rm =T)
sum(dat$ipw_7, na.rm =T)
sum(dat$ipw_8, na.rm =T)
sum(dat$ipw_9, na.rm =T)


# redo weights without colon ----------------------------------------------


# 
# #m1 consists of EN, GI, CO. Avg of 3 for Female, avg of gi co for males
# dat <- dat %>%
#   mutate(m1_ipw = ifelse(gender == "F", (ipw_3+ipw_8)/2, (ipw_8)))
# 
# #m2 consists of P2, NH, OV, P1 
# dat <- dat %>%
#   mutate(m2_ipw = ifelse(gender == "F", (ipw_6+ipw_5 + ipw_2 + ipw_7)/4, (ipw_5 + ipw_2 + ipw_7)/3))
# 
# #m3 consists of LR P2 EN LU NH P1 GI CO
# dat <- dat %>%
#   mutate(m3_ipw = ifelse(gender == "F", (ipw_1+ipw_2+ipw_3+ipw_4+ipw_5+ipw_7+ipw_8)/7,
#                          (ipw_1+ipw_2+ipw_4+ipw_5+ipw_7+ipw_8)/6))
# 
# #m4 consists of LR P2 EN LU OV P1 GI CO
# dat <- dat %>%
#   mutate(m4_ipw = ifelse(gender == "F", (ipw_1+ipw_2+ipw_3+ipw_4+ipw_6+ipw_7+ipw_8)/7,
#                          (ipw_1+ipw_2+ipw_4+ipw_7+ipw_8)/5))
# 
# #m5 LR P2 EN LU NH OV P1 GI CO
# dat <- dat %>%
#   mutate(m5_ipw = ifelse(gender == "F", (ipw_1+ipw_2+ipw_3+ipw_4+ipw_5+ipw_6+ipw_7+ipw_8)/8,
#                          (ipw_1+ipw_2+ipw_4+ipw_5+ipw_7+ipw_8)/6))
# 
# #m6 P1
# dat <- dat %>%
#   mutate(m6_ipw = ipw_2)
# 
# #m7 P2 LU NH OV P1
# dat <- dat %>%
#   mutate(m7_ipw = ifelse(gender == "F", (ipw_2+ipw_5+ipw_6+ipw_7)/4,
#                          (ipw_2+ipw_5+ipw_7)/3))
# 
# #m8 P2 LU NH OV P1 GI
# dat <- dat %>%
#   mutate(m8_ipw = ifelse(gender == "F", (ipw_2+ipw_5+ipw_6+ipw_7+ipw_8)/5,
#                          (ipw_2+ipw_5+ipw_7+ipw_8)/4))
# 
# #m9 LU NH OV P1 GI
# dat <- dat %>%
#   mutate(m9_ipw = ifelse(gender == "F", (ipw_5+ipw_6+ipw_7+ipw_8)/4,
#                          (ipw_5+ipw_7+ipw_8)/3))
# 
# #m10 P2 LU NH P1
# dat <- dat %>%
#   mutate(m10_ipw = (ipw_2+ipw_4+ipw_5+ipw_7)/4)
# 
# 
# #m11 P1 GI
# dat <- dat %>%
#   mutate(m11_ipw = (ipw_7+ipw_9)/2)
# 
# #m12 LU NH OV P1
# dat <- dat %>%
#   mutate(m12_ipw = ifelse(gender == "F", (ipw_4+ipw_5+ipw_6+ipw_7)/4,
#                           (ipw_4+ipw_5+ipw_7)/3))
# 
# #m13 LR EN LU NH OV P1 GI CO
# dat <- dat %>%
#   mutate(m13_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_6+ipw_7+ipw_8)/7,
#                           (ipw_1+ipw_4+ipw_5+ipw_7+ipw_8)/5))
# 
# #m14 GI
# dat <- dat %>%
#   mutate(m14_ipw = ipw_8)
# 
# #m15 LR EN LU NH OV GI CO
# dat <- dat %>%
#   mutate(m15_ipw = ifelse(gender == "F", (ipw_1+ipw_3+ipw_4+ipw_5+ipw_6+ipw_8)/6,
#                           (ipw_1+ipw_4+ipw_5+ipw_8)/4))
# 
# #m16 LU NH OV GI
# dat <- dat %>%
#   mutate(m16_ipw = ifelse(gender == "F", (ipw_4+ipw_5+ipw_6+ipw_8)/4,
#                           (ipw_4+ipw_5+ipw_8)/3))
# 
# #m17 LR P2 EN LU P1 GI CO
# dat <- dat %>%
#   mutate(m17_ipw = ifelse(gender == "F", (ipw_1+ipw_2+ipw_3+ipw_4+ipw_7+ipw_8)/6,
#                           (ipw_1+ipw_2+ipw_4+ipw_7+ipw_8)/5))
# 
# #m18 EN CO
# dat <- dat %>% 
#   mutate(m18_ipw = ifelse(gender == "F", (ipw_3)/1, NA))



# master list -------------------------------------------------------------


#master list of markers x appropriate weights to pick
master <- cbind.data.frame(c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18),
                           c(rep(1, length(m1)),
                             rep(2, length(m2)),
                             rep(3, length(m3)),
                             rep(4, length(m4)),
                             rep(5, length(m5)),
                             rep(6, length(m6)),
                             rep(7, length(m7)),
                             rep(8, length(m8)),
                             rep(9, length(m9)),
                             rep(10, length(m10)),
                             rep(11, length(m11)),
                             rep(12, length(m12)),
                             rep(13, length(m13)),
                             rep(14, length(m14)),
                             rep(15, length(m15)),
                             rep(16, length(m16)),
                             rep(17, length(m17)),
                             rep(18, length(m18))))
names(master) <- c("marker", "ipw_cat")

#add a column that tells if we should run linear or Tobit
#(run Tobit when >10% samples are below LLOD or for SCD30 where 
# 71% are above the ULOD)
master <- master %>% arrange(marker)
master$below <- belowLLOD[belowLLOD$marker != "m_csf_outside_limits", "percent"]
master$model <- ifelse(master$below < .10, 1,
                       ifelse(master$marker == "SCD40L",3,
                              ifelse(master$below < .90, 2, NA)))

master <- master[!is.na(master$model),]
#for the tobit models (model = 2) a censoring value is needed
#find these LLODs and put them in the data frame

#get the minimum observed concentrations for markers
censVal <- dat %>%
  summarise(across(contains("_pg_ml"), .fns = min, na.rm = T))
names(censVal) <- str_remove_all(toupper(names(censVal)), "_PG_ML")

master$min <- (censVal %>% t())

#also add a max col
censVal <- dat %>%
  summarise(across(contains("_pg_ml"), .fns = max, na.rm = T))
names(censVal) <- str_remove_all(toupper(names(censVal)), "_PG_ML")

master$max <- (censVal %>% t())

#actually need to add pg_ml so names match columns in original dat
master$marker <- paste0(tolower(master$marker), "_pg_ml")




# Genotyping --------------------------------------------------------------

#for now, only read in the European since we're limiting analysis to white individuals
# gzList <- list.files("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed", recursive = T)
# gzList <- paste0("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\", gzList)
# gzList <- gzList[str_detect(gzList, "European")]
# #only keep chromosomes that we want ( 1, 4, 7, 9, 12, 17, 18, 19)
# gzList <- gzList[str_detect(gzList, "chr1-")| str_detect(gzList, "chr4-")| str_detect(gzList, "chr7-")|
#                  str_detect(gzList, "chr9-")| str_detect(gzList, "chr12-")| str_detect(gzList, "chr17-")|
#                 str_detect(gzList, "chr18-")| str_detect(gzList, "chr19-")]
# gzList <- gzList[!str_detect(gzList, ".info")]
# gzList <- str_replace_all(gzList, "/", "\\\\")
# 
# gen1 <- map(.x = gzList, .f = possibly(pegas::read.vcf,NA_character_))
# gen2 <- gen1[!is.na(gen1)]
# 
# #function to add the PLCO id onto each data frame in the mapped list
# addPLCOID <- function(x){
#   bind_cols(x, "plco_id" =substr(rownames(x), 1, 10))
# }
# 
# gen2 <- map(gen2, addPLCOID)
# 
# #missing the last 4 chromosomes for people genotyped on the last platform
# #bc of R's memory limits with the pegas::read.vcf function
# #save work, restart to  run these last few
# 
# c19 <- pegas::read.vcf("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\Oncoarray\\European\\chr9-filtered.dose.vcf.gz")
# c7 <- pegas::read.vcf("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\Oncoarray\\European\\chr7-filtered.dose.vcf.gz")
# c4 <- pegas::read.vcf("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\Oncoarray\\European\\chr4-filtered.dose.vcf.gz")
# c9 <- pegas::read.vcf("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\Oncoarray\\European\\chr19-filtered.dose.vcf.gz")
# c18 <- pegas::read.vcf("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\Imputed\\Post_Imputation_QCed\\Oncoarray\\European\\chr18-filtered.dose.vcf.gz")
# 
# 
# gen3 <- map(.x = list(c19,c7,c4,c9, c18), .f = addPLCOID)
# gen2 <- c(gen2,gen3)
# 
# #bind columns for each successive set of 8 elements in gen2 (these are same participants, multiple chromosomes)
# set1 <- bind_cols(gen2[1:7])
# set2 <- bind_cols(gen2[8:14])
# set3 <- bind_cols(gen2[15:21])
# set4 <- bind_cols(gen2[22:28])
# set5 <- bind_cols(gen2[29:36])
# set6 <- bind_cols(gen2[37:44])
# set7 <- bind_cols(gen2[45:52])
# set8 <- bind_cols(gen2[53:60])
# 
# #make ONE dataset by binding rows
# gen.all <- bind_rows(set1, set2, set3, set4, set5, set6, set7, set8)
# gen.all <- gen.all[,!str_detect(names(gen.all), "plco_id")]
# gen.all <- addPLCOID(gen.all)
# rownames(gen.all) <- NULL
# 
# #add an indicator for which platform the participant was genotyped on
# gen.all$platform <- c(rep("GSA",(nrow(set5)+nrow(set4)+nrow(set3)+nrow(set1)+nrow(set2))),
#                       rep("Omni25",nrow(set6)),
#                       rep("OmniX",nrow(set7)),
#                       rep("Oncoarray",nrow(set8)))
# 
# gen.all %>%
#   summarise_all(.funs = ~sum(is.na(.) | (.== ".|.") | (. == ".//.")))

#write.csv(gen.all, "genData.csv")
gen.all <- read.csv("genData.csv")

#merge the genetic data in with the other (mip and dietary and demographics)
#4428 in main
#3986 in genetic data
dat.all <- dat %>%
  inner_join(gen.all, by = "plco_id") #3974
dat <- dat.all

#rename columns in data such that they are the snps rather than chromosome locations
#read in snp list

snplist <- read.delim("snp\\snp.list.updated.txt")

for (i in 1:nrow(snplist)){
  if(str_replace_all(paste0("chr", snplist$Chromosome, ".", snplist$MolecularLocation,
                            str_sub(trimws(snplist$rs_id), -4,-1)), ":", ".")[i] %in% names(dat.all)){
    names(dat.all)[which(str_replace_all(paste0("chr", snplist$Chromosome, ".", snplist$MolecularLocation,
                                                str_sub(trimws(snplist$rs_id), -4,-1)), ":", ".")[i] == names(dat.all))] <- snplist$SNP[i] 
    
  }
}

dat <- dat.all

dat <- dat %>%
  rename(rs12075 = `chr1.159205564.G.A`,
         rs1058396 = `chr18.45739554.G.A`,
         rs8176746 = `chr9.133255935.G.T`,
         rs8176719 = `chr9.133257521.T.TC`)

# for FUT2 (secretor)
dat$rs601338 <- ifelse(is.na(dat$rs601338), 999, dat$rs601338)

#simple dom/recessive
dat <- dat %>%
  mutate(sec = case_when(rs601338 == "A|A" | rs601338 == "A/A" ~ "Nonsecretor",
                         rs601338 %in% c("G/A", "G/G", "A/G", "A|G", "G|A", "G|G") ~ "Secretor",
                         TRUE ~ NA_character_))
#additive model
dat <- dat %>%
  mutate(seh = case_when(sec == "Nonsecretor"~ 0,
                         rs601338 %in% c("G/A", "A/G", "A|G", "G|A") ~ 1,
                         sec == "Secretor" ~ 2,
                         TRUE ~ NA_real_))


#for FUT3 (Lewis) type use method in Cakir 2002
#rs3745635 (T) = 508; rs778986 = 314; rs28362459 =59; rs3894326 = 1067
#rs812936 = 202
# #https://link.springer.com/content/pdf/10.1007/s00277-002-0508-x.pdf


#rs28362459 count number G/C
#rs3894326 count number T-->A
#rs812936 count number G
# dat <- dat %>%
#   mutate(lew = case_when(str_count(rs28362459, "C") == 0 & str_count(rs3894326, "T") == 0 &
#                            str_count(rs812936, "G") == 0 ~ "Wildtype",
#                          str_count(rs3894326, "T") == 2 & str_count(rs812936, "G") == 2 ~ "Mutant",
#                          str_count(rs28362459, "C") == 1 & str_count(rs3894326, "T") == 1 & str_count(rs812936, "G") == 1 ~ "Mutant",
#                          str_count(rs28362459, "C") == 0 & (str_count(rs3894326, "T") == 1 | str_count(rs812936, "G") == 1) ~ "Mutant",
#                          str_count(rs28362459, "C") == 1 & str_count(rs3894326, "T") == 0 & str_count(rs812936, "G") == 1 ~ "Mutant",
#                         TRUE ~ "Subnormal"))


#make lewis alleles by pasting relevant letters
dat$lewisallele1 <- paste0(substr(dat$rs28362459,1,1),
                           substr(dat$rs3894326,1,1),
                           substr(dat$rs3745635,1,1))

dat$lewisallele2 <- paste0(substr(dat$rs28362459,3,3),
                           substr(dat$rs3894326,3,3),
                           substr(dat$rs3745635,3,3))

#see https://www.jbc.org/content/269/46/29271.long
dat$lewisallele1_type <- ifelse(substr(dat$lewisallele1,1,2) == "CT", "le1",
                                ifelse((substr(dat$lewisallele1,1,1)=="C" & substr(dat$lewisallele1,3,3) == "T"), "le2",
                                       ifelse(!substr(dat$lewisallele1,1,1) == ".", "Le", 
                                              ifelse(substr(dat$lewisallele1,2,2) == "T" |
                                                       substr(dat$lewisallele1,3,3) == "T", "le", NA))))
dat$lewisallele2_type <- ifelse(substr(dat$lewisallele2,1,2) == "CT", "le1",
                                ifelse((substr(dat$lewisallele2,1,1)=="C" & substr(dat$lewisallele2,3,3) == "T"), "le2",
                                       ifelse(!substr(dat$lewisallele2,1,1) == ".", "Le", 
                                              ifelse(substr(dat$lewisallele2,2,2) == "T" |
                                                       substr(dat$lewisallele2,3,3) == "T", "le", NA))))

dat$lewis2allele <- paste0(dat$lewisallele1_type, dat$lewisallele2_type)

dat <- dat %>%
  mutate(lew = case_when(str_detect(lewisallele1_type, "L") | str_detect(lewisallele2_type, "L") ~ "Normal",
                         lewis2allele == "NANA" ~ NA_character_,
                         TRUE ~ "Mutant"))

#for additive model
dat <- dat %>%
  mutate(leh = case_when(str_count(lewis2allele, "L")== 2 ~ 2,
                         str_count(lewis2allele, "L") == 1 ~ 1,
                         str_count(lewis2allele, "l") == 2 ~ 0,
                         T ~ NA_real_))


#Duffy Type: rs12075: het = A,B, 2 = A, 0 = B
#rs34599082 = fyx; rs2814778 = null
# #0 duffy null (CC) individuals
dat <- dat %>%
  mutate(duf = case_when(rs34599082 == "C|T" | rs34599082 == "T|C"~"FYA+B-" ,
                         rs2814778 == "C|T" & rs12075 == "A|G" ~ "FYA-B+",
                         rs2814778 == "T|C" & rs12075 == "G|A" ~ "FYA+B-",
                         rs12075 %in% c("A|G", "G|A", "A/G", "G/A") ~ "FYA+B+",
                         rs12075 %in% c("G|G", "G/G") ~ "FYA+B-",
                         rs12075 %in% c("A|A", "A/A") ~ "FYA-B+",
                         TRUE ~ NA_character_))

#Colton Blood Type: rs28362692
dat <- dat %>%
  mutate(col = case_when(rs28362692 %in% c("T|C", "C|T") ~ "CoA+B+",
                         rs28362692 %in% c("C|C", "C/C") ~ "CoA+B-",
                         rs28362692 %in% c("T|T") ~ "CoA-B+",
                         TRUE ~ NA_character_))

#Lutheran type:
dat <- dat %>%
  mutate(lut = case_when(str_count(rs28399653, "G") == 2 ~ "LuA-B+",
                         str_count(rs28399653, "G") == 1 ~ "LuA+B+",
                         TRUE ~ NA_character_))

#Auberger type:
dat <- dat %>%
  mutate(aub = case_when(str_count(rs1135062, "A") == 2 ~ "AuA+B-",
                         str_count(rs1135062, "A") == 1 ~ "AuA+B+",
                         str_count(rs1135062, "G") == 2 ~ "AuA-B+",
                         TRUE ~NA_character_))


#Knops type: (note that rs1047660 and rs1047661 also code Knops but we only 4 and 13 not dominant)
dat <- dat %>%
  mutate(kno = case_when(str_count(rs41274768, "G") == 2 ~ "KnA+B-",
                         str_count(rs41274768, "G") == 1 ~ "KnA+B+",
                         str_count(rs41274768, "A") == 2 ~ "KnA-B+",
                         TRUE ~ NA_character_))

#Scianna (only 6 not dom)

#diego (rs2285644) only 14 not dom; not enough to do anything with

#kidd type
dat <- dat %>%
  mutate(kid = case_when(str_count(rs1058396, "G") == 2 ~ "JkA+B-",
                         str_count(rs1058396, "G") == 1 ~ "JkA+B+",
                         str_count(rs1058396, "A") == 2 ~ "JkA-B+",
                         TRUE ~ NA_character_))


#kell type(s)
dat <- dat %>%
  mutate(kel = case_when(str_count(rs8176058, "G") == 2 ~ "k+k+",
                         str_count(rs8176058, "G") == 1 ~ "K+k+",
                         str_count(rs8176058, "A") == 2 ~ "K+K+",
                         TRUE ~ NA_character_))
#kp antigen (rs8176059) didn't have enough people-> 25 that were kp(a+b+)


#RhE
dat <- dat %>%
  mutate(rhe = case_when(str_count(rs609320, "C") == 2 ~ "EE",
                         str_count(rs609320, "G") == 1 ~ "Ee",
                         str_count(rs609320, "G") == 2 ~ "ee",
                         TRUE ~ NA_character_))

#3 category lewis a/b based on lewis/secretor together
dat <- dat %>%
  mutate(le3 = case_when(lew == "Mutant" ~ "A-B-",
                         sec == "Nonsecretor" ~ "A+B-",
                         sec == "Secretor" ~ "A-B+"))

# #Dombrock type:
#note only 2 hets for rs28362797 and 1 for rs28362798
dat <- dat %>%
  mutate(dom = case_when(str_count(rs11276, "C") == 1 ~ "DoA+B+",
                         str_count(rs11276, "C") == 2 ~ "DoA+B-",
                         str_count(rs11276, "T") == 2 ~ "DoA-B+",
                         T~ NA_character_))


#to make coding ABO allele easier, paste together strings for the two alleles
#gather the needed abo snps
#rs56392308, rs8176747, rs8176746, rs8176743, rs8176741, rs7853989, rs1053878, rs8176720, rs8176719
dat$osnp <- str_replace_all(dat$rs8176719, "TC", "I")
dat$osnp <- str_replace_all(dat$osnp, "T", "D")


dat$aboallele1 <- paste0(substr(dat$osnp,1,1),
                         substr(dat$rs8176720,1,1),
                         substr(dat$rs7853989,1,1),
                         substr(dat$rs1053878,1,1),
                         substr(dat$rs8176741,1,1),
                         substr(dat$rs8176743,1,1),
                         substr(dat$rs8176746,1,1),
                         substr(dat$rs8176747,1,1),
                         substr(dat$rs56392308,1,1))
dat$aboallele2 <- paste0(substr(dat$osnp,3,3),
                         substr(dat$rs8176720,3,3),
                         substr(dat$rs7853989,3,3),
                         substr(dat$rs1053878,3,3),
                         substr(dat$rs8176743,3,1),
                         substr(dat$rs8176746,3,3),
                         substr(dat$rs8176747,3,3),
                         substr(dat$rs56392308,3,3))


dat$aboallele1_type <- ifelse(substr(dat$aboallele1,1,1)=="D" , "O",
                              ifelse(substr(dat$aboallele1,1,1)==".", NA,
                                     ifelse(substr(dat$aboallele1,2,2)=="C" , "B",
                                            ifelse(substr(dat$aboallele1,4,4) == "A", "A2",
                                                   ifelse(!is.na(dat$aboallele1), "A1",NA)))))


dat$aboallele2_type <- ifelse(substr(dat$aboallele2,1,1)=="D" , "O",
                              ifelse(substr(dat$aboallele2,1,1)==".", NA,
                                     ifelse(substr(dat$aboallele2,2,2)=="C", "B",
                                            ifelse(substr(dat$aboallele2,4,4) == "A", "A2",
                                                   ifelse(!is.na(dat$aboallele2), "A1", NA)))))


dat <- dat %>%
  mutate(abo = case_when(aboallele1_type=="O" & aboallele2_type == "O"~"OO",
                         (aboallele1_type=="A1" & aboallele2_type == "A1")~"A1A1",
                         (aboallele1_type=="A2" & aboallele2_type == "A2")~"A2A2",
                         (aboallele1_type=="B" & aboallele2_type == "B")~"BB",
                         (aboallele1_type=="A1" & aboallele2_type == "A2")|(aboallele2_type=="A1" & aboallele1_type == "A2")~"A1A2",
                         (aboallele1_type=="O" & aboallele2_type == "A1")| (aboallele2_type=="O" & aboallele1_type == "A1")~"A1O",
                         (aboallele1_type=="O" & aboallele2_type == "A2")| (aboallele2_type=="O" & aboallele1_type == "A2")~"A2O",
                         (aboallele1_type=="O" & aboallele2_type == "B")| (aboallele2_type=="O" & aboallele1_type == "B")~"BO",
                         (aboallele1_type=="A1" & aboallele2_type == "B")| (aboallele2_type=="A1" & aboallele1_type == "B")~"A1B",
                         (aboallele1_type=="A2" & aboallele2_type == "B")| (aboallele2_type=="A2" & aboallele1_type == "B")~"A2B",
                         TRUE ~ NA_character_))

dat <- dat %>%
  mutate(ab3 = ifelse(abo == "OO", "O",
                      ifelse(abo %in% c("BO", "BB"), "B",
                             ifelse(abo %in% c("A2B", "A1B"), "AB",
                                    ifelse(!is.na(abo), "A", NA)))))

bloodtypelist <- list("ab3", "sec", "lew", "le3", "dom", "rhe", "kel", "kid",
                      "kno", "aub", "lut", "col", "duf")

#recode some of the blood groups with less variation into fewer categories
dat <- dat %>%
  mutate(kno = ifelse(kno == "KnA-B+", "KnA+B+", kno), #5 people
         col = ifelse(col == "CoA-B+", "CoA+B+", col), #3 people
         kel = ifelse(kel == "K+K+", "K+k+", kel)) #4 people


dat$sec <- factor(dat$sec, levels = c("Secretor", "Nonsecretor"))
dat$lew <- factor(dat$lew, levels = c("Normal", "Mutant"))
dat$ab3 <- factor(dat$ab3, levels = c("O", "A", "AB", "B"))
dat$kid <- factor(dat$kid, levels = c("JkA+B+", "JkA+B-", "JkA-B+"))
dat$col <- factor(dat$col, levels = c("CoA+B+", "CoA+B-"))
dat$aub <- factor(dat$aub, levels = c("AuA+B-", "AuA+B+", "AuA-B+"))
dat$le3 <- factor(dat$le3, levels = c("A-B+", "A-B-", "A+B-"))
dat$dom <- factor(dat$dom, levels = c("DoA+B+", "DoA+B-","DoA-B+"))
dat$rhe <- factor(dat$rhe, levels = c("EE", "Ee", "ee"))
dat$abo.all <- dat$abo
dat$abo <- factor(ifelse(dat$abo == "BB", "BO",
                         ifelse(dat$abo == "A2B", "A1B",
                                ifelse(dat$abo == "A2A2", "A2O",
                                       ifelse(dat$abo == "A1A2" | dat$abo == "A1A1", "A1O", dat$abo)))))
dat <- dat %>%
  mutate(a1a = ordered(factor(str_count(paste0(aboallele1_type, aboallele2_type), "A1"))),
         a2a = ordered(factor(str_count(paste0(aboallele1_type, aboallele2_type), "A2"))),
         b1a = ordered(factor(str_count(paste0(aboallele1_type, aboallele2_type), "B"))),
         aaa = ordered(factor(str_count(paste0(aboallele1_type, aboallele2_type), "A"))), #any A
         aba = ordered(factor(2 - str_count(paste0(aboallele1_type, aboallele2_type), "O")))) #any non-O

dat$sec3 <- ordered(factor(dat$sec3))
dat$leh <- ordered(factor(dat$leh))

bloodtypelist <- list("ab3", "sec","le3", "dom", "rhe", "kel", "kid",
                      "kno", "aub", "lut", "col", "duf")


dat <- dat %>%
  mutate(coh = ordered(factor(ifelse(dat$col == "CoA+B+", 1,
                                     ifelse(dat$col == "CoA+B-", 0, NA)))),
         luh = ordered(factor(ifelse(dat$lut == "LuA+B+", 1,
                                     ifelse(dat$lut == "LuA-B+", 0, NA)))),
         knh = ordered(factor(ifelse(dat$kno == "KnA+B+", 1,
                                     ifelse(dat$kno == "KnA+B-", 0, NA)))),
         keh = ordered(factor(ifelse(dat$kel == "K+k+", 1,
                                     ifelse(dat$kel == "k+k+", 0, NA)))),
         doh = ordered(factor(ifelse(dat$dom == "DoA+B+", 1,
                                     ifelse(dat$dom == "DoA+B-", 0,
                                            ifelse(dat$dom == "DoA-B+", 2, NA))))),
         rhh =  ordered(factor(ifelse(dat$rhe == "Ee", 1,
                                      ifelse(dat$rhe == "EE", 0,
                                             ifelse(dat$dom == "ee", 2, NA))))),
         kih = ordered(factor(ifelse(dat$kid == "JkA+B+", 1,
                                     ifelse(dat$kid == "JkA+B-", 0,
                                            ifelse(dat$kid == "JkA-B+", 2, NA))))),
         auh = ordered(factor(ifelse(dat$aub == "AuA+B+", 1,
                                     ifelse(dat$aub == "AuA+B-", 0,
                                            ifelse(dat$aub == "AuA-B+", 2, NA))))),
         duh = ordered(factor(ifelse(dat$duf == "FYA+B+", 1,
                                     ifelse(dat$duf == "FYA+B-", 0,
                                            ifelse(dat$duf == "FYA-B+", 2, NA)))))
  )

#dat$abo <- factor(dat$abo, levels = c("OO", "BO", "A1O", "A2O", "A1B"),ordered = F)
dat$duf <- factor(dat$duf, levels = c("FYA+B-", "FYA+B+", "FYA-B+"), ordered = F)

dat$leb <- factor(ifelse(dat$le3 == "A+B-", "LeA+",
                         ifelse(dat$le3 == "A-B+", "B+", NA)))

dat$ab3 <- factor(dat$ab3, levels = c("O", "A", "AB", "B"))

dat$sec <- relevel(factor(dat$sec, ordered = F), ref = "Secretor")
dat$dom <- factor(dat$dom, levels = c("DoA+B-", "DoA+B+", "DoA-B+"))
dat$rhe <- factor(dat$rhe, levels = c("EE", "Ee", "ee"))
dat$kid <- factor(dat$kid, levels = c("JkA+B-", "JkA+B+", "JkA-B+"))
dat$aub <- factor(dat$aub, levels = c("AuA+B-", "AuA-B+", "AuA+B+"))

#dominant list:
domlist <- list("sec", "dom", "rhe", "kid", "aub", "duf", "leb", "ab3")

#additive list: le2, seh, aba, aaa, a1a, a2a, b1a, coh, knh, keh, doh, rhh, auh, luh, duh
addlist <- list('le2', 'seh', 'a1a', 'a2a', 'b1a', 'aaa', 'aba', 'kih',
                'coh', 'knh', 'keh', 'doh', 'rhh', 'auh', 'luh', 'duh')




# Read in Principal Components --------------------------------------------
# pc <- read.delim("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\PrincipalComponents\\European\\GSA.step7.evec",
#                  sep = "\t")
# pc$id <- trimws(pc$id)
# 
# pc$id <- str_split(pc$id, "    ")
# for (i in 2:11){
#   for (j in 1:nrow(pc)){
#     pc[j,i] <-  pc$id[[j]][i]
#   }
# }
# 
# pc2 <- read.delim("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\PrincipalComponents\\European\\Omni5.step7.evec",
#                  sep = "\t")
# pc2$id <- trimws(pc2$id)
# 
# pc2$id <- str_split(pc2$id, "    ")
# for (i in 2:11){
#   for (j in 1:nrow(pc2)){
#     pc2[j,i] <-  pc2$id[[j]][i]
#   }
# }
# 
# pc3 <- read.delim("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\PrincipalComponents\\European\\Omni25.step7.evec",
#                   sep = "\t")
# pc3$id <- trimws(pc3$id)
# 
# pc3$id <- str_split(pc3$id, "    ")
# for (i in 2:11){
#   for (j in 1:nrow(pc3)){
#     pc3[j,i] <-  pc3$id[[j]][i]
#   }
# }
# 
# pc4 <- read.delim("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\PrincipalComponents\\European\\OmniX.step7.evec",
#                   sep = "\t")
# pc4$id <- trimws(pc4$id)
# 
# pc4$id <- str_split(pc4$id, "    ")
# for (i in 2:11){
#   for (j in 1:nrow(pc4)){
#     pc4[j,i] <-  pc4$id[[j]][i]
#   }
# }
# 
# pc5 <- read.delim("C:\\Users\\vanalstensc\\Documents\\Rcode\\MIP_BloodType\\PrincipalComponents\\European\\Oncoarray.step7.evec",
#                   sep = "\t")
# pc5$id <- trimws(pc5$id)
# 
# pc5$id <- str_split(pc5$id, "    ")
# for (i in 2:11){
#   for (j in 1:nrow(pc5)){
#     pc5[j,i] <-  pc5$id[[j]][i]
#   }
# }
# 
# names(pc) <- c("plco_id", paste0("pc",1:10))
# names(pc2) <- names(pc3) <- names(pc4) <- names(pc5) <- names(pc)
# 
# pc$array <- "GSA"
# pc2$array <- "Omni5"
# pc3$array <- "Omni25"
# pc4$array <- "OmniX"
# pc5$array <- "Oncoarray"
# 
# pc <- bind_rows(pc, pc2, pc3, pc4, pc5)
# 
# for (i in 1:nrow(pc)){
#   pc$plco_id[i] <- pc$plco_id[[i]][1]
# }
# 
# #make sure PCs are numeric
# pc <- pc %>%
#   mutate(across(2:11, .fns = ~(as.numeric(as.character(trimws(.))))))
# 
# class(pc$plco_id) <- "character"
# 
# length(unique(pc$plco_id))
# write.csv(pc, "pc1.csv)

#takes forever- write out to file so can read this in diretly
pc <- read_csv("pc1.csv")

#for sanity incase the join doesn't work
#dat.copy <- dat

dat <- dat %>%
  left_join(pc, by = "plco_id")


# Binary Logistic Function ---------------------------------------------------------
#Makers with >50% missing just dichotomize
overHalfMiss <- belowLLOD[belowLLOD$percent>.5 & belowLLOD$percent<.9,"marker"]
overHalfMiss <- str_replace_all(overHalfMiss, "_outside_limits", "_pg_ml")


#binary need to be 0/1
dat <- dat %>%
  mutate(across(ends_with("_bin"), .fns = ~(.-1)))



dat2 <- dat %>%
  select_at(overHalfMiss) %>%
  mutate(across(.fns = ~ifelse(. == min(., na.rm = T), 1, 2)))
names(dat2) <- str_replace_all(names(dat2), "_pg_ml", "_bin")

#SCD40L also bc 70% above ULOD
dat3 <- dat %>%
  select(scd40l_pg_ml) %>%
  mutate(across(.fns = ~(ifelse(. >=3.9, 2,1))))
names(dat3) <- str_replace_all(names(dat3), "_pg_ml", "_bin")

dat <- bind_cols(dat, dat2, dat3)


#names of markers with 1/2 missing
halfMiss <- dat %>%
  select_at(vars(ends_with("_bin"))) %>%
  names()

dat$overall_study <- factor(dat$overall_study)


#lot number to exclude depending on marker/panel from CRC study
master$exclude <- 0
master$exclude <- ifelse(master$marker %in% c("adipsin_pg_ml_pg_ml", "sap_pg_ml_pg_ml", "crp_pg_ml_pg_ml"), 2395931,
                         ifelse(master$marker %in% c("segfr_pg_ml_pg_ml", "sgp_130_pg_ml_pg_ml", "sil4r_pg_ml_pg_ml", "sil6r_pg_ml_pg_ml",
                                                     "silrii_pg_ml_pg_ml", "stnfri_pg_ml_pg_ml", "stnfrii_pg_ml_pg_ml", "svegfr2_pg_ml_pg_ml",
                                                     "svegfr3_pg_ml_pg_ml"), 2343687,
                                ifelse(master$marker %in% c("egf_pg_ml_pg_ml", "eotaxin_pg_ml_pg_ml", "fgf_2_pg_ml_pg_ml", "g_csf_pg_ml_pg_ml",
                                                            "gro_pg_ml_pg_ml", "vegf_pg_ml_pg_ml", "tnfa_pg_ml_pg_ml", "tnf_b_pg_ml_pg_ml",
                                                            "tgf_a_pg_ml_pg_ml", "ifng_pg_ml_pg_ml", "ifna2_pg_ml_pg_ml", "il_10_pg_ml_pg_ml",
                                                            "il_1ra_pg_ml_pg_ml", "il_7_pg_ml_pg_ml", "il_8_pg_ml_pg_ml", "ip_10_pg_ml_pg_ml", "il_33_pg_ml_pg_ml",
                                                            "il_29_ifnl1_pg_ml_pg_ml"), 2343697,
                                       ifelse(master$marker %in% c("resistin_pg_ml_pg_ml", "mcp_1_pg_ml_pg_ml"), 2069501,
                                              ifelse(master$marker %in% c("x6ckine_pg_ml_pg_ml", "bca_1_pg_ml_pg_ml", "ctack_pg_ml_pg_ml",
                                                                          "ena_78_pg_ml_pg_ml", "eotaxin_2_pg_ml_pg_ml", "eotaxin_3_pg_ml_pg_ml", "mcp_2_pg_ml_pg_ml",
                                                                          "mcp_4_pg_ml_pg_ml", "mip_1d_pg_ml_pg_ml", "scf_pg_ml_pg_ml", "sdf_1a_b_pg_ml_pg_ml", "trail_pg_ml_pg_ml",
                                                                          "tarc_pg_ml_pg_ml", "tslp_pg_ml_pg_ml", "tpo_pg_ml_pg_ml"), 2344094,
                                                     ifelse(master$marker %in% c("ccl20_mip3a_pg_ml_pg_ml", "ccl19_mip3b_pg_ml_pg_ml", "cxcl9_mig_pg_ml_pg_ml",
                                                                                 "cxcl11_i_tac_pg_ml_pg_ml", "cxcl6_gcp2_pg_ml_pg_ml"), 2344152, 0))))))

#write.csv(dat, "stabilized_ipw_0911.csv")
#function to run logistic regression
runLogReg <- function(marker, bloodtype){
  
  #first: extract the appropriate weights to apply
  ipw_cat <- paste0("m", master[master$marker == str_replace(marker, "_bin", "_pg_ml_pg_ml"), "ipw_cat"], "_ipw")
  ipw_cat <- dat[, names(dat) == ipw_cat]
  
  #are we excluding any crc?
  ex <- master[master$marker == str_replace(marker, "_bin", "_pg_ml_pg_ml"), "exclude"]
  
  #subset to covariates and the relevant marker and blood type
  #note I'm not including overall_study due to linear dependence issues
  
  df <- dat[, names(dat) %in% c(marker, bloodtype, "gender", "overall_study",
                                "age5", "cig_stat", "pc1", "pc2", "pc3",
                                "pc4", "pc5", "C6P3_LotNum", "AD5P1_LotNum", "C16P2_LotNum",
                                "CVD3P2_LotNum", "SR9P1_LotNum", "C17P1_LotNum")]
  
  df$ipw_category <- ipw_cat
  df$ipw_category <- as.numeric(df$ipw_category)
  
  #depending on marker, exclude some of the people in CRC based on batch
  if (ex == 2344152){
    df <- df[is.na(df$C6P3_LotNum) | !(df$C6P3_LotNum == 2344152),]
  } else if (ex == 2344094){
    df <- df[is.na(df$C16P2_LotNum) | !(df$C16P2_LotNum == 2344094),]
  } else if (ex == 2069501){
    df <- df[is.na(df$AD5P1_LotNum) | !(df$AD5P1_LotNum == 2069501),]
  } else if (ex == 2343697){
    df <- df[is.na(df$C17P1_LotNum) | !(df$C17P1_LotNum == 2343697),]
  } else if (ex == 2343687){
    df <- df[is.na(df$SR9P1_LotNum) | !(df$SR9P1_LotNum == 2343687),]
  } else if (ex == 2395931){
    df <- df[is.na(df$CVD3P2_LotNum) | !(df$CVD3P2_LotNum == 2343687),]
  }
  
  #model formula
  fmla <- formula(paste0(marker, "~", bloodtype, "+ gender + age5 + cig_stat",
                         "+ pc1 + pc2 + pc3 + pc4 + pc5 + overall_study"))
  
  df <- df[!is.na(df$ipw_category),]
  
  #run the model
  mod <- glm(fmla, data = df, family = binomial(link = "logit")) 
  mod.sum <- summary(mod)
  wald <- survey::regTermTest(mod, bloodtype, method = "Wald")$p
  
  #tidy it up in a nice data frame
  lmt <- lmtest::coeftest(mod, vcov = vcov(mod))
  
  return(list(lmt, wald, mod))
}


#function to extract relevant results from the logistic models
getLogRes <- function(mod, wald, mod2, problem){
  
  mark <- str_split(mod2$terms[1], " ")[[2]]
  
  res <- cbind.data.frame(mod[2:5,],
                          rep(wald, 4),
                          rep(mark, 4),
                          rep(problem, 4))
}



# Run Logistic Regressions ------------------------------------------------
#binary need to be 0/1
dat <- dat %>%
  mutate(across(ends_with("_bin"), .fns = ~(.-1)))

#DOMINANT GENETIC MODELS
cr <- cross_df(list(arg1 = halfMiss, arg2 = domlist))

cr <- cr[-210,]#removes the one that results in an error
#logregs <- map2(.x = cr$arg1, .y = cr$arg2, .f = possibly(runLogReg, NA))
logregs <- map2(.x = cr$arg1, .y = cr$arg2, .f = quietly(.f = runLogReg))
#logregs <- map2(.x = cr$arg1, .y = cr$arg2, .f = runLogReg)

# find the elements which produced warnings
problem <- which(str_detect(logregs %>% map("warnings"), "fitted"))

#extract models
logisticmodlist <- logregs %>% map("result")

#extract the wald tests
waldlist <- sapply(logisticmodlist, "[", 2)
#extract coeftest
modlist <- sapply(logisticmodlist, "[", 1)
#extract the actual model objects
modlist2 <- sapply(logisticmodlist, "[", 3)

#new problem vector
prob2 <- c(1:255)
prob2 <- ifelse(prob2 %in% problem, "Did Not Converge", "Converged")

logregdf <- pmap_dfr(.l = list(modlist, waldlist, modlist2, prob2), .f = possibly(getLogRes, tibble()))
logregdf$term <- rownames(logregdf)

logregdf <- logregdf %>% filter(!str_detect(term, "gender") &
                                  !str_detect(term, "cig_stat") & !str_detect(term, "pc") &
                                  !str_detect(term, "age5"))



# Joint Additive for ABO --------------------------------------------------


#for ABO, also do a set of models where #a1, a2, b alleles all in model together
runABOLog <- function(marker){
  
  #first: extract the appropriate weights to apply
  ipw_cat <- paste0("m", master[master$marker == str_replace(marker, "_bin", "_pg_ml_pg_ml"), "ipw_cat"], "_ipw")
  ipw_cat <- dat[, names(dat) == ipw_cat]
  
  #are we excluding any crc?
  ex <- master[master$marker == str_replace(marker, "_bin", "_pg_ml_pg_ml"), "exclude"]
  
  #subset to covariates and the relevant marker and blood type
  #note I'm not including overall_study due to linear dependence issues
  
  df <- dat[, names(dat) %in% c(marker, "a1a","a2a","b1a", "gender", "overall_study",
                                "age5", "cig_stat", "pc1", "pc2", "pc3",
                                "pc4", "pc5", "C6P3_LotNum", "AD5P1_LotNum", "C16P2_LotNum",
                                "CVD3P2_LotNum", "SR9P1_LotNum", "C17P1_LotNum")]
  
  df$ipw_category <- ipw_cat
  df$ipw_category <- as.numeric(df$ipw_category)
  
  #model formula
  fmla <- formula(paste0(marker, "~ a1a + a2a + b1a + gender + age5 + cig_stat",
                         "+ pc1 + pc2 + pc3 + pc4 + pc5 + overall_study"))
  
  #depending on marker, exclude some of the people in CRC based on batch
  if (ex == 2344152){
    df <- df[is.na(df$C6P3_LotNum) | !(df$C6P3_LotNum == 2344152),]
  } else if (ex == 2344094){
    df <- df[is.na(df$C16P2_LotNum) | !(df$C16P2_LotNum == 2344094),]
  } else if (ex == 2069501){
    df <- df[is.na(df$AD5P1_LotNum) | !(df$AD5P1_LotNum == 2069501),]
  } else if (ex == 2343697){
    df <- df[is.na(df$C17P1_LotNum) | !(df$C17P1_LotNum == 2343697),]
  } else if (ex == 2343687){
    df <- df[is.na(df$SR9P1_LotNum) | !(df$SR9P1_LotNum == 2343687),]
  } else if (ex == 2395931){
    df <- df[is.na(df$CVD3P2_LotNum) | !(df$CVD3P2_LotNum == 2343687),]
  }
  
  
  #run the model
  mod <- glm(fmla, data = df, weights = ipw_category, family = binomial(link = "logit"))
  mod.sum <- summary(mod)
  wald <- survey::regTermTest(mod, ~a1a+a2a+b1a, method = "Wald")$p
  
  #tidy it up in a nice data frame
  lmt <- lmtest::coeftest(mod, vcov = vcov(mod))
  
  return(list(lmt, wald, mod))
  
}


abolog <- map(.x = halfMiss, .f = quietly(runABOLog))
# find the elements which produced warnings
problem <- which(str_detect(abolog %>% map("warnings"), "fitted"))
#extract models
logisticmodlist <- abolog %>% map("result")
#extract the wald tests
waldlist <- sapply(logisticmodlist, "[", 2)
#extract coeftest
modlist <- sapply(logisticmodlist, "[", 1)
#extract the actual model objects
modlist2 <- sapply(logisticmodlist, "[", 3)

#new problem vector
prob2 <- c(1:32)
prob2 <- ifelse(prob2 %in% problem, "Did Not Converge", "Converged")

abologdf <- pmap_dfr(.l = list(modlist, waldlist, modlist2, prob2), .f = possibly(getLogRes, tibble()))
abologdf$term <- rownames(abologdf)

abologdf <- abologdf %>% filter(!str_detect(term, "gender") &
                                  !str_detect(term, "cig_stat") & !str_detect(term, "pc") &
                                  !str_detect(term, "age5"))



# Linear Regression Function ------------------------------------------------------

#markers not in the halfMiss list
underHalfMiss <- belowLLOD[belowLLOD$percent<.5 & belowLLOD$marker != "scd40l_outside_limits","marker"]
underHalfMiss <- underHalfMiss[underHalfMiss != "m_csf_pg_ml"]
underHalfMiss <- str_replace_all(underHalfMiss, "_outside_limits", "_pg_ml")

#function to run logistic regression
runLinReg <- function(marker, bloodtype){
  
  #first: extract the appropriate weights to apply
  ipw_cat <- paste0("m", master[master$marker == str_replace(marker, "_pg_ml", "_pg_ml_pg_ml"), "ipw_cat"], "_ipw")
  ipw_cat <- dat[, names(dat) == ipw_cat]
  
  
  
  #are we excluding any crc?
  ex <- master[master$marker == str_replace(marker, "_pg_ml", "_pg_ml_pg_ml"), "exclude"]
  
  #subset to covariates and the relevant marker and blood type
  #note I'm not including overall_study due to linear dependence issues
  
  df <- dat[, names(dat) %in% c(marker, bloodtype, "gender", "overall_study",
                                "age5", "cig_stat", "pc1", "pc2", "pc3",
                                "pc4", "pc5", "C6P3_LotNum", "AD5P1_LotNum", "C16P2_LotNum",
                                "CVD3P2_LotNum", "SR9P1_LotNum", "C17P1_LotNum")]
  df$ipw_category <- ipw_cat
  df$ipw_category <- as.numeric(df$ipw_category)
  
  #model formula
  fmla <- formula(paste0(marker, "~", bloodtype, "+ gender + age5 + cig_stat",
                         "+ pc1 + pc2 + pc3 + pc4 + pc5 + overall_study"))
  
  #depending on marker, exclude some of the people in CRC based on batch
  if (ex == 2344152){
    df <- df[is.na(df$C6P3_LotNum) | !(df$C6P3_LotNum == 2344152),]
  } else if (ex == 2344094){
    df <- df[is.na(df$C16P2_LotNum) | !(df$C16P2_LotNum == 2344094),]
  } else if (ex == 2069501){
    df <- df[is.na(df$AD5P1_LotNum) | !(df$AD5P1_LotNum == 2069501),]
  } else if (ex == 2343697){
    df <- df[is.na(df$C17P1_LotNum) | !(df$C17P1_LotNum == 2343697),]
  } else if (ex == 2343687){
    df <- df[is.na(df$SR9P1_LotNum) | !(df$SR9P1_LotNum == 2343687),]
  } else if (ex == 2395931){
    df <- df[is.na(df$CVD3P2_LotNum) | !(df$CVD3P2_LotNum == 2343687),]
  }
  
  df <- df[!is.na(df$ipw_category),]
  
  #run the model
  mod <- lm(fmla, data = df, weights = ipw_category)
  mod.sum <- summary(mod)
  wald <- survey::regTermTest(mod, bloodtype, method = "Wald")$p
  
  #tidy it up in a nice data frame
  lmt <- lmtest::coeftest(mod, vcov = vcov(mod))
  
  return(list(lmt, wald, mod))
}

# #function to extract relevant results from the logistic models
# getLinRes <- function(mod, wald){
#   
#   
#   mod.sum <- broom::tidy(mod, conf.int = TRUE)
#   
#   mod.sum$marker <- names(attr(mod$terms,which = "dataClasses"))[1]
#   mod.sum$waldp <- wald
#   if(nrow(mod.sum) > 23){
#     return(mod.sum[2:10,])
#   } else if(nrow(mod.sum) == 23){
#     return(mod.sum[2:5,])
#   } else if(nrow(mod.sum) == 22){
#     return(mod.sum[2:4,])
#   } else if(nrow(mod.sum) == 21){
#     return(mod.sum[2:3,])
#   } else{
#     return(mod.sum[2:5,])
#   }
#   #result has 15, 13, or 14 rows depending on number of blood type levels
# }

# Run Linear Regs ---------------------------------------------------------

#Dominant Models
cr2 <- cross_df(list(arg1 = underHalfMiss, arg2 = domlist))
abolog <- map2(.x = cr2$arg1, .y = cr2$arg2, .f = possibly(quietly(runLinReg), NA))
abolog <- abolog[!is.na(abolog)]
# find the elements which produced warnings
problem <- which(str_detect(abolog %>% map("warnings"), "fitted"))
#extract models
logisticmodlist <- abolog %>% map("result")
#extract the wald tests
waldlist <- sapply(logisticmodlist, "[", 2)
#extract coeftest
modlist <- sapply(logisticmodlist, "[", 1)
#extract the actual model objects
modlist2 <- sapply(logisticmodlist, "[", 3)

#new problem vector
prob2 <- c(1:456)
prob2 <- ifelse(prob2 %in% problem, "Did Not Converge", "Converged")

linregdf <- pmap_dfr(.l = list(modlist, waldlist, modlist2, prob2), .f = possibly(getLogRes, tibble()))
linregdf$term <- rownames(linregdf)

linregdf <- linregdf %>% filter(!str_detect(term, "gender") &
                                  !str_detect(term, "cig_stat") & !str_detect(term, "pc") &
                                  !str_detect(term, "age5"))


#Additive Models
cr3 <- cross_df(list(arg1 = underHalfMiss, arg2 = addlist))
abolog <- map2(.x = cr3$arg1, .y = cr3$arg2, .f = possibly(quietly(runLinReg), NA))
abolog <- abolog[!is.na(abolog)]
# find the elements which produced warnings
problem <- which(str_detect(abolog %>% map("warnings"), "fitted"))
#extract models
logisticmodlist <- abolog %>% map("result")
#extract the wald tests
waldlist <- sapply(logisticmodlist, "[", 2)
#extract coeftest
modlist <- sapply(logisticmodlist, "[", 1)
#extract the actual model objects
modlist2 <- sapply(logisticmodlist, "[", 3)

#new problem vector
prob2 <- c(1:855)
prob2 <- ifelse(prob2 %in% problem, "Did Not Converge", "Converged")

linregdf2 <- pmap_dfr(.l = list(modlist, waldlist, modlist2, prob2), .f = possibly(getLogRes, tibble()))
linregdf2$term <- rownames(linregdf2)

linregdf2 <- linregdf2 %>% filter(!str_detect(term, "gender") &
                                    !str_detect(term, "cig_stat") & !str_detect(term, "pc") &
                                    !str_detect(term, "age5"))

#Joint additive for abo

#for ABO, also do a set of models where #a1, a2, b alleles all in model together
runABOLin <- function(marker){
  
  #first: extract the appropriate weights to apply
  ipw_cat <- paste0("m", master[master$marker == str_replace(marker, "_pg_ml", "_pg_ml_pg_ml"), "ipw_cat"], "_ipw")
  ipw_cat <- dat[, names(dat) == ipw_cat]
  
  
  
  #are we excluding any crc?
  ex <- master[master$marker == str_replace(marker, "_pg_ml", "_pg_ml_pg_ml"), "exclude"]
  
  #subset to covariates and the relevant marker and blood type
  #note I'm not including overall_study due to linear dependence issues
  
  df <- dat[, names(dat) %in% c(marker,"a1a","a2a","b1a", "gender", "overall_study",
                                "age5", "cig_stat", "pc1", "pc2", "pc3",
                                "pc4", "pc5", "C6P3_LotNum", "AD5P1_LotNum", "C16P2_LotNum",
                                "CVD3P2_LotNum", "SR9P1_LotNum", "C17P1_LotNum")]
  df$ipw_category <- ipw_cat
  df$ipw_category <- as.numeric(df$ipw_category)
  
  #model formula
  fmla <- formula(paste0(marker, "~", "a1a + a2a + b1a + gender + age5 + cig_stat",
                         "+ pc1 + pc2 + pc3 + pc4 + pc5 + overall_study"))
  
  #depending on marker, exclude some of the people in CRC based on batch
  if (ex == 2344152){
    df <- df[is.na(df$C6P3_LotNum) | !(df$C6P3_LotNum == 2344152),]
  } else if (ex == 2344094){
    df <- df[is.na(df$C16P2_LotNum) | !(df$C16P2_LotNum == 2344094),]
  } else if (ex == 2069501){
    df <- df[is.na(df$AD5P1_LotNum) | !(df$AD5P1_LotNum == 2069501),]
  } else if (ex == 2343697){
    df <- df[is.na(df$C17P1_LotNum) | !(df$C17P1_LotNum == 2343697),]
  } else if (ex == 2343687){
    df <- df[is.na(df$SR9P1_LotNum) | !(df$SR9P1_LotNum == 2343687),]
  } else if (ex == 2395931){
    df <- df[is.na(df$CVD3P2_LotNum) | !(df$CVD3P2_LotNum == 2343687),]
  }
  
  
  #model formula
  fmla <- formula(paste0(marker, "~ a1a + a2a + b1a + gender + age5 + cig_stat",
                         "+ pc1 + pc2 + pc3 + pc4 + pc5 + overall_study"))
  
  #run the model
  mod <- lm(fmla, data = df, weights = ipw_category)
  #mod <- lm(fmla, data = df) 
  wald <- survey::regTermTest(mod, ~a1a+a2a+b1a, method = "Wald")$p
  #tidy it up in a nice data frame
  lmt <- lmtest::coeftest(mod, vcov = vcov(mod))
  
  return(list(lmt, wald, mod))
  
}

#function to extract relevant results from the logistic models
getLogRes2 <- function(mod, wald, mod2, problem){
  
  mark <- str_split(mod2$terms[1], " ")[[2]]
  
  res <- cbind.data.frame(mod[2:7,],
                          rep(wald, 6),
                          rep(mark, 6),
                          rep(problem, 6))
}


abolog <- map(.x = underHalfMiss, .f =possibly(quietly(runABOLin), NA))
abolog <- abolog[!is.na(abolog)]
# find the elements which produced warnings
problem <- which(str_detect(abolog %>% map("warnings"), "fitted"))
#extract models
logisticmodlist <- abolog %>% map("result")
#extract the wald tests
waldlist <- sapply(logisticmodlist, "[", 2)
#extract coeftest
modlist <- sapply(logisticmodlist, "[", 1)
#extract the actual model objects
modlist2 <- sapply(logisticmodlist, "[", 3)

#new problem vector
prob2 <- c(1:57)
prob2 <- ifelse(prob2 %in% problem, "Did Not Converge", "Converged")

abolindf <- pmap_dfr(.l = list(modlist, waldlist, modlist2, prob2), .f = possibly(getLogRes2, tibble()))
abolindf$term <- rownames(abolindf)

abolindf <- abolindf %>% filter(!str_detect(term, "gender") &
                                  !str_detect(term, "cig_stat") & !str_detect(term, "pc") &
                                  !str_detect(term, "age5"))


# Summarize ---------------------------------------------------------------
#put the linear reg + logistic reg results together for dominant
linlog.dom <- bind_rows(linregdf, logregdf)
linlog.dom$logistic <- c(rep(F, nrow(linregdf)), rep(T, nrow(logregdf)))
linlog.dom$bloodtypecat <- substr(linlog.dom$term, 1,3)
linlog.dom$waldp <- linlog.dom$`rep(wald, 4)`

linlog.dom.simp <- linlog.dom %>%
  filter(!bloodtypecat == "abo")

12*64


linregdf %>%
  mutate(waldp = `rep(wald, 4)`)%>%
  mutate(bonferroni = ifelse(waldp < 0.05/768, 1, 0)) %>%
  filter(bonferroni==1)

linregdf2 %>%
  mutate(waldp = `rep(wald, 4)`)%>%
  mutate(bonferroni = ifelse(waldp < 0.05/768, 1, 0)) %>%
  filter(bonferroni==1)


length(unique(linlog.dom.simp$waldp))
#bonferroni p-value to compare to
linlog.dom <- linlog.dom %>%
  mutate(bonferroni = ifelse(waldp < 0.05/680, T, F))

length(unique(linlog.dom$waldp))
#bonferroni p-value to compare to
linlog.dom.simp <- linlog.dom.simp %>%
  mutate(bonferroni = ifelse(waldp < 0.05/680, T, F))

table(linlog.dom$bonferroni)

signif.dom <- linlog.dom.simp[linlog.dom.simp$bonferroni == T,]
#write.csv(linlog.dom, "dominant09102020.csv")

#write.csv(linlog.dom, "dominant08272020.csv")

#QQ type plots separated out for each blood type
gg_qqplot <- function(ps, bloodtype, ci = 0.95) {
  ps <- unique(ps)
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 19, size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())+
    ggtitle(label = bloodtype)
}


ptests <- linlog.dom.simp %>% filter(bloodtypecat=="ab3") %>% dplyr::select(waldp, bloodtypecat)
abp <- gg_qqplot(ps = ptests$waldp, "ABO")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="duf") %>% dplyr::select(waldp)
dup <- gg_qqplot(ps = ptests$waldp, "Duffy")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="rhe") %>% dplyr::select(waldp)
rhp <- gg_qqplot(ps = ptests$waldp, "RhE")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="aub") %>% dplyr::select(waldp)
aup <- gg_qqplot(ps = ptests$waldp, "Auberger")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="kid") %>% dplyr::select(waldp)
kip <- gg_qqplot(ps = ptests$waldp, "Kidd")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="sec") %>% dplyr::select(waldp)
sep <- gg_qqplot(ps = ptests$waldp, "Secretor")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="dom") %>% dplyr::select(waldp)
dop <- gg_qqplot(ps = ptests$waldp, "Dombrock")
ptests <- linlog.dom.simp %>% filter(bloodtypecat=="leb") %>% dplyr::select(waldp)
lep <- gg_qqplot(ps = ptests$waldp, "Lewis")

library(gridExtra)

gA <- ggplotGrob(abp)
gB <- ggplotGrob(dup)
gD <- ggplotGrob(rhp)
gC <- ggplotGrob(aup)
gE <- ggplotGrob(kip)
gF <- ggplotGrob(sep)
gG <- ggplotGrob(dop)
gH <- ggplotGrob(lep)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5],
                           gC$widths[2:5], gD$widths[2:5],
                           gE$widths[2:5], gF$widths[2:5],
                           gG$widths[2:5], gH$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
gD$widths[2:5] <- as.list(maxWidth)
gE$widths[2:5] <- as.list(maxWidth)
gF$widths[2:5] <- as.list(maxWidth)
gG$widths[2:5] <- as.list(maxWidth)
gH$widths[2:5] <- as.list(maxWidth)

library(patchwork)

tiff("SupplementaryFigure1.pdf")
(abp + dup + rhp)/(aup + kip + sep)/ (dop + lep + plot_spacer()) +
  plot_annotation(tag_levels = "A", title = "Supplementary Figure 1.",
                  caption = "QQ plots of blood type-marker associations for dominant models.\n95% bootstrapped confidence interval is shown in the dotted line,\nand only points > -log10(7.5) were statistically significant with\nBonferroni correction.")
dev.off()

pdf("dom_pp_plots_09102020.pdf")
grid.arrange(grobs = list(gA, gB, gC, gD, gE, gF, gG, gH),
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6),
                                   c(7,8,NA)),
             widths = c(2,2,2))
dev.off()


#put the linear reg + logistic reg results together for additive

linlog.add <- bind_rows(linregdf2,
                        logregdf,
                        abolindf, abologdf)

linlog.add$logistic <- c(rep(F, nrow(linregdf2)), 
                         rep(T, nrow(logregdf)),
                         rep(F, nrow(abolindf)), rep(T, nrow(abologdf)))


#adjust the p-values (is bonferroni sig?)
linlog.add <- linlog.add %>%
  filter(term != "a1a.Q" & term != " a2a.Q" & term != "aaa.L" & term != "aba.L" & term != "b1a.Q" &
           term != "a2a.Q" & term != "genderM" & term != "cig_stat" & term != "age5" & term != "pc1" & term != "aaa.Q" &
           term != "aba.Q" & term != "auh.Q" & term != "doh.Q" & term != "duh.Q" & term != "kih.Q")
linlog.add <- linlog.add %>%
  filter(!str_detect(term, "age") &!str_detect(term, "gender") & !str_detect(term, "cig_stat") )

linlog.add <- linlog.add %>%
  filter(!str_detect(term , ".Q"))

linlog.add$bloodtypecat <- substr(linlog.add$term, 1,3)

linlog.add <- linlog.add %>%
  filter(!bloodtypecat == "ab3")

linlog.add$waldp <- linlog.add$`rep(wald, 4)`
length(unique(linlog.add$waldp)) 
table(linlog.add$waldp == 0)


linlog.add <- linlog.add %>%
  mutate(bonferroni = ifelse(waldp < 0.05/(1079) & !is.na(waldp), T, F))

table(linlog.add$bonferroni)
signif.add <- linlog.add[linlog.add$bonferroni,]

#write.csv(linlog.add, "additive09102020.csv")
#write.csv(linlog.add, "additive09052020.csv")
#write.csv(linlog.add, "additive08272020.csv")
#write.csv(linlog.dom, "dominant08272020.csv")

#make more QQ-ish p plots
ptests <- linlog.add %>% filter(bloodtypecat=="a1a") %>% dplyr::select(waldp, bloodtypecat)
a1p <- gg_qqplot(ps = ptests$waldp, "ABO:A1")
ptests <- linlog.add %>% filter(bloodtypecat=="a2a") %>% dplyr::select(waldp, bloodtypecat)
a2p <- gg_qqplot(ps = ptests$waldp, "ABO:A2")
ptests <- linlog.add %>% filter(bloodtypecat=="b1a") %>% dplyr::select(waldp, bloodtypecat)
b1p <- gg_qqplot(ps = ptests$waldp, "ABO:B1")
ptests <- linlog.add %>% filter(bloodtypecat=="duh") %>% dplyr::select(waldp)
dup <- gg_qqplot(ps = ptests$waldp, "Duffy")
ptests <- linlog.add %>% filter(bloodtypecat=="rhh") %>% dplyr::select(waldp)
rhp <- gg_qqplot(ps = ptests$waldp, "RhE")
ptests <- linlog.add %>% filter(bloodtypecat=="auh") %>% dplyr::select(waldp)
aup <- gg_qqplot(ps = ptests$waldp, "Auberger")
ptests <- linlog.add %>% filter(bloodtypecat=="kih") %>% dplyr::select(waldp)
kip <- gg_qqplot(ps = ptests$waldp, "Kidd")
ptests <- linlog.add %>% filter(bloodtypecat=="seh") %>% dplyr::select(waldp)
sep <- gg_qqplot(ps = ptests$waldp, "Secretor")
ptests <- linlog.add %>% filter(bloodtypecat=="doh") %>% dplyr::select(waldp)
dop <- gg_qqplot(ps = ptests$waldp, "Dombrock")
ptests <- linlog.add %>% filter(bloodtypecat=="le2") %>% dplyr::select(waldp)
lep <- gg_qqplot(ps = ptests$waldp, "Lewis")


gA <- ggplotGrob(a1p)
gA1 <- ggplotGrob(a2p)
gAb <- ggplotGrob(b1p)
gB <- ggplotGrob(dup)
gD <- ggplotGrob(rhp)
gC <- ggplotGrob(aup)
gE <- ggplotGrob(kip)
gF <- ggplotGrob(sep)
gG <- ggplotGrob(dop)
gH <- ggplotGrob(lep)
maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5], gA1$widths[2:5], gAb$widths[2:5],
                           gC$widths[2:5], gD$widths[2:5],
                           gE$widths[2:5], gF$widths[2:5],
                           gG$widths[2:5], gH$widths[2:5])
gA$widths[2:5] <- as.list(maxWidth)
gAb$widths[2:5] <- as.list(maxWidth)
gA1$widths[2:5] <- as.list(maxWidth)
gB$widths[2:5] <- as.list(maxWidth)
gC$widths[2:5] <- as.list(maxWidth)
gD$widths[2:5] <- as.list(maxWidth)
gE$widths[2:5] <- as.list(maxWidth)
gF$widths[2:5] <- as.list(maxWidth)
gG$widths[2:5] <- as.list(maxWidth)
gH$widths[2:5] <- as.list(maxWidth)

pdf("add_pp_plots.pdf")
grid.arrange(grobs = list(gA, gA1, gAb, gB, gC, gD, gE, gF, gG, gH),
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6),
                                   c(7,8,9),
                                   c(10,11,NA)),
             widths = c(2,2,2))
dev.off()

library(patchwork)

pdf("additive_pp_significant.pdf")
(a1p + a2p)/ (b1p + dup)
dev.off()

# #write these results out
# write.csv(dat, "analysis_data.csv")
# write.csv(linlog.dom, "dominant.csv")
# write.csv(linlog.add, "additive.csv")


table(cor(dat[, names(dat) %in% signif.dom$marker],
          use = "pairwise.complete")<.3)

# Assumption Checks -------------------------------------------------------

#function to do some assumption checking: residual normality
residNorm <- function(mod.resid){
  
  shap <- shapiro.test(mod.resid$.resid)
  return(shap$p.value)
}

residHist <- function(mod.resid){
  mod.resid %>%
    ggplot()+
    geom_histogram(aes(x=.resid))
}


#function to do some assumption checking: homoscedasticity
homosced <- function(mod.resid){
  
  names(mod.resid)[2] <- "marker"
  # mod.resid %>%
  #   ggplot() +
  #   geom_point(aes(x = .resid, y = .fitted))
  return(lmtest::bptest(mod.resid$marker~mod.resid$.resid)$p.value)
}





augs <- map(linmodlist, .f = possibly(broom::augment, NA_real_))
norms <- map(augs, possibly(residNorm, NA_real_))
norms <- unlist(norms)

r <- map(augs, residHist)
hom <- map(augs, possibly(homosced, NA_real_))

pdf("residualsHist.pdf")
r
dev.off()

table(unlist(hom) < 0.05)
# Supplementary Tables 3+4 ------------------------------------------------
linlog.dom <- read.csv("dominant09052020.csv")
linlog.add <- read.csv("additive09052020.csv")

sup3 <- 
  linlog.dom %>%
  mutate(`Estimate (95% CI)` = paste0(sprintf(estimate, fmt ="%.2f"), " (",
                                      sprintf(conf.low, fmt ="%.2f"), " - ",
                                      sprintf(conf.high, fmt ="%.2f"), ")")) %>%
  mutate(`Blood Group` = case_when(bloodtypecat == "leb" ~ "Lewis",
                                   bloodtypecat %in% c("abo", "aba", "a1a", "b1a", "a2a") ~ "ABO",
                                   bloodtypecat == "rhe" ~ "RhE",
                                   bloodtypecat == "aub" ~ "Aub",
                                   bloodtypecat == "kel" ~ "Kell",
                                   bloodtypecat == "kid" ~ "Kidd",
                                   bloodtypecat == "duf" ~ "Duffy",
                                   bloodtypecat == "sec" ~ "Secretor",
                                   bloodtypecat == "lut" ~ "Lutheran",
                                   bloodtypecat == "kno" ~ "Knops",
                                   bloodtypecat == "dom" ~ "Dombrock",
                                   bloodtypecat == "col" ~ "Colton")) %>%
  mutate(`Antigen Phenotype` = str_remove(term, bloodtypecat)) %>%
  mutate(`Antigen Phenotype` = ifelse(`Antigen Phenotype` == "BO", "B",
                                      ifelse(`Antigen Phenotype` == "A1O", "A1",
                                             ifelse(`Antigen Phenotype` == "A2O", "A2",
                                                    ifelse(`Antigen Phenotype` == "A1B", "AB", `Antigen Phenotype`))))) %>%  
  mutate(Marker = toupper(str_replace_all(str_remove_all(str_remove_all(marker, "_pg_ml"), "_bin"), "_", " "))) %>%
  mutate(`P-value` = sprintf(p.value, fmt = "%.3f")) %>%
  mutate(`Wald Joint P-value` = sprintf(waldp, fmt = "%.3f")) %>%
  mutate(`Bonferroni Significant?` = bonferroni) %>%
  relocate(1, .after = 4) %>%
  mutate(`Estimate (95% CI)` = ifelse(str_detect(`Estimate (95% CI)`, "NA"), "--",
                                      `Estimate (95% CI)`)) %>%
  mutate(`Estimate (95% CI)` = ifelse(nchar(`Estimate (95% CI)`)>25, "--",
                                      `Estimate (95% CI)`)) %>%
  mutate(`Estimate (95% CI)` = ifelse(str_detect(`Estimate (95% CI)`, "Inf"), "--",
                                      `Estimate (95% CI)`))


its.lit <- 
  linlog.add %>%
  mutate(estimate = ifelse(is.infinite(conf.high) | is.na(conf.high), NA, estimate),
         conf.low = ifelse(is.infinite(conf.high) | is.na(conf.high), NA, conf.low),
         conf.high = ifelse(is.infinite(conf.high) | is.na(conf.high), NA, conf.high)) %>%
  mutate(`Estimate (95% CI)` = paste0(sprintf(estimate, fmt ="%.2f"), " (",
                                      sprintf(conf.low, fmt ="%.2f"), " - ",
                                      sprintf(conf.high, fmt ="%.2f"), ")")) %>%
  mutate(`Blood Group` = case_when(bloodtypecat == "le2" ~ "Lewis",
                                   bloodtypecat %in% c("aaa", "aba", "a1a", "b1a", "a2a") ~ "ABO",
                                   bloodtypecat == "rhh" ~ "RhE",
                                   bloodtypecat == "auh" ~ "Aub",
                                   bloodtypecat == "keh" ~ "Kell",
                                   bloodtypecat == "kih" ~ "Kidd",
                                   bloodtypecat == "duh" ~ "Duffy",
                                   bloodtypecat == "seh" ~ "Secretor",
                                   bloodtypecat == "luh" ~ "Lutheran",
                                   bloodtypecat == "knh" ~ "Knops",
                                   bloodtypecat == "doh" ~ "Dombrock",
                                   bloodtypecat == "coh" ~ "Colton")) %>%
  mutate(`Allele` = case_when(bloodtypecat == "le2" ~ "Null",
                              bloodtypecat %in% c("aaa") ~ "Any A",
                              bloodtypecat %in% c("aba") ~ "Any Non-O",
                              bloodtypecat %in% c( "a1a") ~ "A1",
                              bloodtypecat %in% c("a2a") ~ "A2",
                              bloodtypecat %in% c("b1a") ~ "B",
                              bloodtypecat == "rhh" ~ "E",
                              bloodtypecat == "auh" ~ "B",
                              bloodtypecat == "keh" ~ "K",
                              bloodtypecat == "kih" ~ "B",
                              bloodtypecat == "duh" ~ "B",
                              bloodtypecat == "seh" ~ "Secretor",
                              bloodtypecat == "luh" ~ "B",
                              bloodtypecat == "knh" ~ "B",
                              bloodtypecat == "doh" ~ "B",
                              bloodtypecat == "coh" ~ "B")) %>%
  mutate(Marker = toupper(str_replace_all(str_remove_all(str_remove_all(marker, "_pg_ml"), "_bin"), "_", " "))) %>%
  mutate(`P-value` = sprintf(p.value, fmt = "%.3f")) %>%
  mutate(`Wald Joint P-value` = sprintf(waldp, fmt = "%.3f")) %>%
  mutate(`Bonferroni Significant?` = bonferroni) %>%
  mutate(`Estimate (95% CI)` = ifelse(str_detect(`Estimate (95% CI)`, "NA"), "--",
                                      `Estimate (95% CI)`))
# Write Supplemental Material to an Excel File ----------------------------

#make data frame of study characteristics
studyChar <- rbind.data.frame(cases = c(526, 526, 301, 149, 171, 284, 63),
                              controls = c(592, 625, 301, 149, 344, 284, 63),
                              percent_baseline = c(100,100,100, 11.4,87,90.1,80.2),
                              median_time_to_dx = c("2.9 (1.1 - 5.1)", "3.7 (2.5 - 5.9)",
                                                    "8.0 (5.0 - 13.9)", "4.2 (2.8 - 6.7)",
                                                    "Not Reported", "5.3 (2.1 - 9.1)", "6.5 (3.6 - 9.5)"),
                              inclusion_criteria_1 = rep("Blood Serum Available",7),
                              inclusion_criteria_2 = rep("Consent", 7),
                              inclusion_criteria_3 = rep("No Cancer History", 7),
                              inclusion_crieria_4 = c("Complete Smoking History", "Complete Smoking History", "",
                                                      "Female", "Complete Smoking History", "Female", ""),
                              inclusion_criteria_5 = c("No multiple cancers during follow up",
                                                       "No multiple cancers during follow up", "", 
                                                       "No control with oopherectomy", "Non-Hispanic White",
                                                       "No control with hysterectomy", ""),
                              inclusion_criteria_6 = c(rep("", 4), "Received Sigmoidoscopy", "", ""),
                              inclusion_criteria_7 = c(rep("", 4), "No self reported ulcerative colitis, Crohn's Gardner's or familial polyposis",
                                                       "", ""),
                              matching_criteria_1 = c(rep("5 year age category", 4), "Used Lung, NHL, Ovarian samples",
                                                      "5 year age category", "1 year age"),
                              matching_criteria_2 = c("Sex", "Sex", "Race", "Race", "", "Race", "Sex"),
                              matching_criteria_3 = c("Randomization Year", "Study Year of Blood draw",
                                                      "Study Center", "Study Center", "", "Study Year of Blood draw",
                                                      "Race"),
                              matching_criteria_4 = c("Smoking Status", "Smoking Status", "Time of blood draw (am/pm)",
                                                      "Time of blood draw (am/pm)", "", "Randomization Year", "Study Year of Exit"),
                              matching_criteria_5 = c("Cumulative Smoking at baseline", "Cumulative Smoking at baseline",
                                                      "Date", "Date", "", "", "Number Freeze/thaw cycles"),
                              matching_criteria_6 = c("Time since quitting", "Time since quitting", rep("",5)))
names(studyChar) <- c("Lung", "LungRep", "NHL", "Ovarian", "Colorectal", "Endometrial", "UpperGI")


#output out supplemental material for formatting
writexl::write_xlsx(x = list("concentrationsMeans" = mean_sd, "concentrationQuartiles" = markerConc,
                             "partial_correlations" = as.data.frame(partial.corr),
                             "study_characteristics" = studyChar, "SuppTable3" = sup3,
                             "SuppTable4" = its.lit),
                    path = "mip_output_nhw_09062020.xlsx", col_names = TRUE, format_headers = TRUE)


writexl::write_xlsx(its.lit, "supp4_09052020.xlsx")
writexl::write_xlsx(sup3, "supp3_09052020.xlsx")

# Barplots Showing Concentrations in Different Blood Groups ---------------

makeLineRangePlot <- function(mark, bt){
  
  d<- dat %>%
    dplyr::select(duf, ab3, svegfr2_pg_ml, svegfr3_pg_ml, sgp130_pg_ml,
                  mcp_1_pg_ml, mcp_4_pg_ml, gro_pg_ml, eotaxin_pg_ml,
                  ena_78_pg_ml, tarc_pg_ml, cxcl6_gcp2_pg_ml, duh, b1a, a1a) %>%
    select_at(c(mark, bt))
  
  names(d) <- c("mark", "bt")
  
  log10Pg <- expression(paste("Log"[10], plain( pg/ml)))
  
  title <- trimws(toupper(str_remove_all(str_replace_all(mark, "_", " "), "pg ml")))
  
  title <- ifelse(title == "CXCL6 GCP2", "CXCL6/GCP2",
                  ifelse(title == "ENA 78", "CXCL5/ENA78",
                         ifelse(title == "EOTAXIN", "CCL11/EOTAXIN",
                                ifelse(title == "GRO", "CXCL1/GRO",
                                       ifelse(title == "MCP 1", "CCL2/MCP1",
                                              ifelse(title == "MCP 4", "CCL13/MCP4",
                                                     ifelse(title == "TARC", "CCL17/TARC", title)))))))
  
  d %>%
    drop_na() %>%
    ggplot(aes(x = bt, y= mark)) +
    stat_summary(fun.data = "mean_cl_normal", color = ifelse(bt == "duf" | bt == "duh", "red", "blue"),
                 size = .2, geom= "pointrange") +
    labs(x = "", y = log10Pg) +
    ggtitle(title) + theme(panel.grid.minor = element_blank())+
    theme_bw()# +
  # theme(axis.text.x = element_text(face = "bold", size = 12),
  #       axis.text.y = element_text(face = "plain", size = 12))
}



p1 <- makeLineRangePlot("tarc_pg_ml", "duf")
p2 <- makeLineRangePlot("mcp_1_pg_ml", "duf")
p3 <- makeLineRangePlot("mcp_4_pg_ml", "duf")
p4 <- makeLineRangePlot("gro_pg_ml", "duf")
p5 <- makeLineRangePlot("cxcl6_gcp2_pg_ml", "duf")
p6 <- makeLineRangePlot("ena_78_pg_ml", "duf")
p7 <- makeLineRangePlot("eotaxin_pg_ml", "duf")
p8 <- makeLineRangePlot("svegfr3_pg_ml", "ab3")
p9 <- makeLineRangePlot("svegfr2_pg_ml", "ab3")
p10 <- makeLineRangePlot("sgp130_pg_ml", "ab3")


(p1 + p2 + p3)/(p4 + p5 + p6)/(p7 + plot_spacer() + plot_spacer())



dat %>%
  ggplot(aes(x = duf, y = tarc_pg_ml))+ geom_boxplot()


#plots for additive models
signif.add <- signif.add[6:nrow(signif.add),]

p1 <- makeLineRangePlot("tarc_pg_ml", "duh")
p2 <- makeLineRangePlot("mcp_1_pg_ml", "duh")
p3 <- makeLineRangePlot("mcp_4_pg_ml", "duh")
p4 <- makeLineRangePlot("gro_pg_ml", "duh")
p5 <- makeLineRangePlot("cxcl6_gcp2_pg_ml", "duh")
p6 <- makeLineRangePlot("ena_78_pg_ml", "duh")
p7 <- makeLineRangePlot("eotaxin_pg_ml", "duh")
p8 <- makeLineRangePlot("svegfr3_pg_ml", "ab3")
p9 <- makeLineRangePlot("svegfr2_pg_ml", "ab3")
p10 <- makeLineRangePlot("sgp130_pg_ml", "ab3")




# Estimated Marginal Means ------------------------------------------------

library(emmeans)
m1 <- lm(tarc_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)
em1 <- emmeans(m1, specs = "duf")
sem1 <- summary(em1)

m2 <- lm(mcp_1_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em2 <- emmeans(m2, specs = "duf")
sem2 <- summary(em2)


m3 <- lm(mcp_4_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em3 <- emmeans(m3, specs = "duf")
sem3 <- summary(em3)


m4 <- lm(gro_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em4<- emmeans(m4, specs = "duf")
sem4 <- summary(em4)

m5 <- lm(cxcl6_gcp2_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em5 <- emmeans(m5, specs = "duf")
sem5 <- summary(em5)


m6 <- lm(ena_78_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em6<- emmeans(m6, specs = "duf")
sem6 <- summary(em6)


m7 <- lm(eotaxin_pg_ml ~ duf  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em7 <- emmeans(m7, specs = "duf")
sem7 <- summary(em7)


m8 <- lm(svegfr3_pg_ml ~ ab3  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em8 <- emmeans(m8, specs = "ab3")
sem8 <- summary(em8)

m9 <- lm(svegfr2_pg_ml ~ ab3  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = dat)

em9 <- emmeans(m9, specs = "ab3")
sem9 <- summary(em9)

m10 <- lm(sgp130_pg_ml ~ ab3  + gender + age5 + cig_stat+ pc1 +
            pc2 + pc3 + pc4 + pc5 + overall_study,
          data = dat)

em10 <- emmeans(m10, specs = "ab3")
sem10 <- summary(em10)


ggeffects::ggeffect(m10, terms = "abonew")

str(em10)
#data frame with raw and estimated marginal means that can be plotted
duf.sum <- dat %>%
  group_by(duf) %>%
  drop_na(duf) %>%
  mutate(n = n())%>%
  group_by(duf)%>%
  summarise(meanTarc = mean(tarc_pg_ml, na.rm =T),
            meanMCP1 = mean(mcp_1_pg_ml, na.rm =T),
            meanMCP4 = mean(mcp_4_pg_ml, na.rm =T),
            meanEotaxin = mean(eotaxin_pg_ml, na.rm =T),
            meanEna78 = mean(ena_78_pg_ml, na.rm=T),
            meanGro = mean(gro_pg_ml, na.rm =T),
            meanCXCL6 = mean(cxcl6_gcp2_pg_ml, na.rm =T),
            seTarc = sd(tarc_pg_ml, na.rm =T)/sqrt(n),
            seMCP1 = sd(mcp_1_pg_ml, na.rm =T)/sqrt(n),
            seMCP4 = sd(mcp_4_pg_ml, na.rm =T)/sqrt(n),
            seEotaxin = sd(eotaxin_pg_ml, na.rm =T)/sqrt(n),
            seEna78 = sd(ena_78_pg_ml, na.rm=T)/sqrt(n),
            seGro = sd(gro_pg_ml, na.rm =T)/sqrt(n),
            seCXCL6 = sd(cxcl6_gcp2_pg_ml, na.rm =T)/sqrt(n),
            .groups = 'keep') 

duf.sum <- duf.sum[c(1,1834,3885),]

abo.sum <- dat %>%
  mutate(abonew = ifelse(aboallele1_type< aboallele2_type, 
                         paste0(aboallele1_type, aboallele2_type),
                         paste0(aboallele2_type, aboallele1_type))) %>%
  mutate(abonew = ifelse(abonew == "BO" | abonew == "BB", "B",
                         ifelse(abonew == "A1O" | abonew == "A1A1", "A1",
                                ifelse(abonew == "A2O" | abonew == "A2A2", "A2",
                                       ifelse(abonew == "OO", "O", abonew)))))# %>%
# group_by(abonew) %>%
# summarise(meanSVEGFR2 = mean(svegfr2_pg_ml, na.rm = T),
#           meanSVEGFR3 = mean(svegfr3_pg_ml, na.rm =T),
#           meanSGP130 = mean(sgp130_pg_ml, na.rm =T),
#           .groups = "keep")

m8 <- lm(svegfr3_pg_ml ~ abonew  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = abo.sum)

em8 <- emmeans(m8, specs = "abonew")
sem8 <- summary(em8)

m9 <- lm(svegfr2_pg_ml ~ abonew  + gender + age5 + cig_stat+ pc1 +
           pc2 + pc3 + pc4 + pc5 + overall_study,
         data = abo.sum)

em9 <- emmeans(m9, specs = "abonew")
sem9 <- summary(em9)

m10 <- lm(sgp130_pg_ml ~ abonew  + gender + age5 + cig_stat+ pc1 +
            pc2 + pc3 + pc4 + pc5 + overall_study,
          data = abo.sum)

em10 <- emmeans(m10, specs = "abonew")
sem10 <- summary(em10)

abo.sum <- abo.sum %>%
  group_by(abonew) %>%
  summarise(meanSVEGFR2 = mean(svegfr2_pg_ml, na.rm = T),
            meanSVEGFR3 = mean(svegfr3_pg_ml, na.rm =T),
            meanSGP130 = mean(sgp130_pg_ml, na.rm =T),
            .groups = "keep")


# Supplementary Figures ---------------------------------------------------
linlog.dom <- read.csv("dominant09052020.csv")
linlog.add <- read.csv("additive09052020.csv")
#QQ type plots separated out for each blood type
gg_qqplot <- function(ps, bloodtype, bonf, bound) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    bonf = bonf,
    clower   = bound[,1],
    cupper   = bound[,2]
  )
  
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed, color = bonf), shape = 19, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = "dashed") +
    geom_ribbon(aes(ymin = clower, ymax = cupper, x = expected),
                alpha = .2, fill = "pink")+
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank())+
    scale_color_manual(values = c("black", "red"), labels = c("No", "Yes"))+
    labs(color = "Significant with\nBonferroni Correction?")+
    scale_x_continuous(breaks =c(0.5,1,1.5,2,2.5))+
    theme(plot.title = element_text(size = 25, face = "bold"),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))+
    ggtitle(label = bloodtype)
}
#function to run simulation
simFun <- function(n = 85){
  x <- runif(n = n, min = 0, max = 1) #p values are uniformly distributed
  x <- sort(x) #sort them in order so when we bind each row = one rank
  x <- tibble(x) #return as tibble
}

set.seed(8675309)
#######################################################################################################
#Auburger
ptests <- linlog.dom %>% filter(bloodtypecat=="aub") %>% dplyr::select(waldp, bonferroni)
ptests <- get_dupes(ptests) %>% group_by(waldp)%>% mutate(count = lag(dupe_count)) %>% filter(is.na(count)) #88
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
aup <- gg_qqplot(ptests$waldp, "Auberger", bonf = ptests$bonferroni, bound = boundsAub)
###########################################################################################################
#Dombrock
ptests <- linlog.dom %>% filter(bloodtypecat=="dom") %>% dplyr::select(waldp, bonferroni)
ptests <- get_dupes(ptests) %>% group_by(waldp)%>% mutate(count = lag(dupe_count)) %>% filter(is.na(count)) #88
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
dop <- gg_qqplot(ptests$waldp, "Dombrock", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################
#Kidd
ptests <- linlog.dom %>% filter(bloodtypecat=="kid") %>% dplyr::select(waldp, bonferroni)
ptests <- get_dupes(ptests) %>% group_by(waldp)%>% mutate(count = lag(dupe_count)) %>% filter(is.na(count)) #88
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
kip <- gg_qqplot(ptests$waldp, "Kidd", bonf = ptests$bonferroni, bound = boundsAub)
###############################################################################################################
#lewis
ptests <- linlog.dom %>% filter(bloodtypecat=="leb") %>% dplyr::select(waldp, bonferroni)
runit <- rep(86, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
lep <- gg_qqplot(ptests$waldp, "Lewis", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################
#rhe
ptests <- linlog.dom %>% filter(bloodtypecat=="rhe") %>% dplyr::select(waldp, bonferroni)
ptests <- get_dupes(ptests) %>% group_by(waldp)%>% mutate(count = lag(dupe_count)) %>% filter(is.na(count)) #88
runit <- rep(88, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
rhp <- gg_qqplot(ptests$waldp, "Rhe", bonf = ptests$bonferroni, bound = boundsAub)
##################################################################################################################
#sec
ptests <- linlog.dom %>% filter(bloodtypecat=="sec") %>% dplyr::select(waldp, bonferroni)
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
sep <- gg_qqplot(ptests$waldp, "Secretor", bonf = ptests$bonferroni, bound = boundsAub)
####################################################################################################################
#luh
ptests <- linlog.add %>% filter(bloodtypecat=="luh") %>% dplyr::select(waldp, bonferroni)
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
lup <- gg_qqplot(ptests$waldp, "Lutheran", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################
#knh
ptests <- linlog.add %>% filter(bloodtypecat=="knh") %>% dplyr::select(waldp, bonferroni)
runit <- rep(86, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
knp <- gg_qqplot(ptests$waldp, "Knops", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################
#keh
ptests <- linlog.add %>% filter(bloodtypecat=="keh") %>% dplyr::select(waldp, bonferroni)
runit <- rep(87, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
kep <- gg_qqplot(ptests$waldp, "Kell", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################
#coh
ptests <- linlog.add %>% filter(bloodtypecat=="coh") %>% dplyr::select(waldp, bonferroni)
runit <- rep(88, 10000) #87 p values generated
testSimAub <- purrr::map_dfc(.x = runit,
                             .f = simFun)
testSimAub <- t(testSimAub)
#get 95% CI of each col
boundsAub <- matrixStats::colQuantiles(testSimAub, probs = c(.025, .975))
boundsAub[,1] <- -log10(boundsAub[,1])
boundsAub[,2] <- -log10(boundsAub[,2])
cop <- gg_qqplot(ptests$waldp, "Colton", bonf = ptests$bonferroni, bound = boundsAub)
################################################################################################################

pdf("supplementary_figures_1_qqplots09052020.pdf")
cop
kep
knp
sep
lep
kip
rhp
aup
dop
dev.off()

#double check that last 5 PCs aren't correlated with the markers
#randomly PC8 was correlated with MCP4 but none of the others were
#significantly correlated
cor.test(dat$pc6, dat$svegfr3_pg_ml)
cor.test(dat$pc7, dat$svegfr3_pg_ml)
cor.test(dat$pc8, dat$svegfr3_pg_ml)
cor.test(dat$pc9, dat$svegfr3_pg_ml)
cor.test(dat$pc10, dat$svegfr3_pg_ml)

cor.test(dat$pc6, dat$svegfr2_pg_ml)
cor.test(dat$pc7, dat$svegfr3_pg_ml)
cor.test(dat$pc8, dat$svegfr2_pg_ml)
cor.test(dat$pc9, dat$svegfr3_pg_ml)
cor.test(dat$pc10, dat$svegfr2_pg_ml)

cor.test(dat$pc6, dat$sgp130_pg_ml)
cor.test(dat$pc7, dat$sgp130_pg_ml)
cor.test(dat$pc8, dat$sgp130_pg_ml)
cor.test(dat$pc9, dat$sgp130_pg_ml)
cor.test(dat$pc10, dat$sgp130_pg_ml)

cor.test(dat$pc6, dat$tarc_pg_ml)
cor.test(dat$pc7, dat$tarc_pg_ml)
cor.test(dat$pc8, dat$tarc_pg_ml)
cor.test(dat$pc9, dat$tarc_pg_ml)
cor.test(dat$pc10, dat$tarc_pg_ml)

cor.test(dat$pc6, dat$eotaxin_pg_ml)
cor.test(dat$pc7, dat$eotaxin_pg_ml)
cor.test(dat$pc8, dat$eotaxin_pg_ml)
cor.test(dat$pc9, dat$eotaxin_pg_ml)
cor.test(dat$pc10, dat$eotaxin_pg_ml)

cor.test(dat$pc6, dat$mcp_4_pg_ml)
cor.test(dat$pc7, dat$mcp_4_pg_ml)
cor.test(dat$pc8, dat$mcp_4_pg_ml)
cor.test(dat$pc9, dat$mcp_4_pg_ml)
cor.test(dat$pc10, dat$mcp_4_pg_ml)

cor.test(dat$pc6, dat$mcp_1_pg_ml)
cor.test(dat$pc7, dat$mcp_1_pg_ml)
cor.test(dat$pc8, dat$mcp_1_pg_ml)
cor.test(dat$pc9, dat$mcp_1_pg_ml)
cor.test(dat$pc10, dat$mcp_1_pg_ml)

cor.test(dat$pc6, dat$cxcl6_gcp2_pg_ml)
cor.test(dat$pc7, dat$cxcl6_gcp2_pg_ml)
cor.test(dat$pc8, dat$cxcl6_gcp2_pg_ml)
cor.test(dat$pc9, dat$cxcl6_gcp2_pg_ml)
cor.test(dat$pc10, dat$cxcl6_gcp2_pg_ml)

cor.test(dat$pc6, dat$ena_78_pg_ml)
cor.test(dat$pc7, dat$ena_78_pg_ml)
cor.test(dat$pc8, dat$ena_78_pg_ml)
cor.test(dat$pc9, dat$ena_78_pg_ml)
cor.test(dat$pc10, dat$ena_78_pg_ml)
cor.test(dat$pc5, dat$ena_78_pg_ml)

cor.test(dat$pc6, dat$gro_pg_ml)
cor.test(dat$pc7, dat$gro_pg_ml)
cor.test(dat$pc8, dat$gro_pg_ml)
cor.test(dat$pc9, dat$gro_pg_ml)
cor.test(dat$pc10, dat$gro_pg_ml)
cor.test(dat$pc5, dat$gro_pg_ml)

#######################################
dclus1<-survey::svydesign(id=~1, weights=~m13_ipw, data=dat[!is.na(dat$m13_ipw),])
summary(dclus1)

(tbl <- svytable(~sch.wide+stype, dclus1))




sum(dat$m13_ipw)
sum(dat$m5_ipw)
sum(dat$m10_ipw)
