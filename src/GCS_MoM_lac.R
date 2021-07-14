# SET ENVIRONMENT
# install.packages(c('remotes', 'here', 'tidyverse', 'RcppArmadillo', 'svglite'))
# remotes::install_github(c('julou/ggCustomTJ', 'hadley/multidplyr'))
# remotes::install_github('vanNimwegenLab/vngMoM', auth_token='xxx')
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cowplot)
  library(RcppArmadillo)
  library(vngMoM)
  library(ggCustomTJ)
})

theme_set(theme_cowplot() + theme(title = element_text(size = rel(1/1.14)),
                                  strip.text.x=element_text(margin=margin(t=1, b=2)),
                                  strip.text.y=element_text(margin=margin(l=2, r=1))) )
theme_cowplot_legend_inset <- function(.rel=0.7) theme(legend.title=element_text(size=rel(.rel), face='bold'),
                                                       legend.text=element_text(size=rel(.rel)))

# set a parallel environment to run multidplyr (ALL packages explicitely loaded before will be loaded too)
library(multidplyr)
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  new_cluster(.options = callr::r_session_options(user_profile = FALSE)) %>% # disable user profile to avoid endless startup with `renv`
  cluster_library( # load currently loaded packages on each core (but multidplyr)
    sessionInfo()$otherPkgs %>% names() %>% setdiff("multidplyr")
  )
# Packages loaded from now on won't be attached to the multidplyr cluster

dir.create(here("slogs"), showWarnings=FALSE) # create a directory to store logs from the queue


# SET VARIABLES
dl <- 0.065                # pixel size (µm)
vertical_cutoff <- 1 / dl  # after it touched this coordinate a cell is discarded
min_growth_rate <- 4e-5    # growth rate threshold to discard non growing cells before the switch (sec-1)
use_kaiser2018_params <- TRUE

data2preproc_dir <- function(.d)
  str_match(.d, '20\\d{6}') %>% na.omit %>% as.character %>% 
  file.path('.', 'preproc', .)
data2preproc_file <- function(.f)
  basename(.f) %>% sub("ExportedCellStats_", "", .) %>% 
  tools::file_path_sans_ext() %>% paste0("_frames.txt")
data2preproc <- function(.f)
  file.path(data2preproc_dir(.f), data2preproc_file(.f))

# EXPERIMENTAL CONDITIONS AND DATA PATHS
myconditions <- list(
  # # CONSTANT CONDITIONS
  # list(condition='mg1655', duration=c(720, 720), dt=180, medium=c('glucose', 'lactose'),
  #      paths=c("./data_MoM_ms/MG1655_glu_lac")),
  # list(condition='glucose', duration=1560, dt=180, medium='glucose',
  #      paths=c("./data_MoM_ms/glucose")),
  list(condition='lactose', duration=1560, dt=180, medium='lactose',
       paths=c("./data_Kaiser2018/lactose")),

  # WITH/WITHOUT GROWTH ARREST
  list(condition='switch_lactulose_TMG20_stdIllum', duration=c(480, 720), dt=360, 
       medium=c('glycerol', 'lactulose+TMG'),
       paths=c(#"./data_thomas/20180704/20180704_glyc_lactulose_TMG20uM_curated/", # with dt=180 slow growth and limited induction
         "data_thomas/20180709/20180709_glyc_lactuloseTMG20uM_curated/", "./data_thomas/20180711/20180711_glyc_lactuloseTMG20uM_curated/")),
  list(condition='switch_lactulose_stdIllum', duration=c(480, 720), dt=360, 
       medium=c('glycerol', 'lactulose'),
       paths=c("./data_thomas/20180710/20180710_glyc_lactulose_curated/")),
  list(condition='switch_lactulose_TMG20', duration=c(480, 720, 480), dt=360, 
       medium=c('glycerol', 'lactulose+TMG', 'glycerol'),
       paths=c("./data_thomas/20181008/20181008_glyc_lactuloseTMG20uM_curated/", "./data_thomas/20181009/20181009_glyc_lactuloseTMG20uM_curated/")),
  list(condition='switch_lactulose', duration=c(480, 720, 480), dt=360, 
       medium=c('glycerol', 'lactulose', 'glycerol'),
       paths=c("./data_thomas/20181024/20181024_glyc_lactulose_curated/")),
  list(condition='switch_glycerol_TMG20', duration=c(480, 720), dt=360, 
       medium=c('glycerol', 'glycerol+TMG'),
       paths=c("./data_thomas/20180712/20180712_glyc_glycTMG20uM_curated/")),
  
  # SUGAR MIXTURES
  list(condition='M9zero', duration=48*60, dt=180, medium='M9zero',
       paths=c("./data_theo/20210708/20210708_curated/20210708_S5")),

  list(condition='glu2uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S0/", "./data_theo/20210305/20210305_curated/20210305_S0/", 
               "./data_theo/20210513/20210513_curated/20210513_S0/")),
  list(condition='glu4uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210513/20210513_curated/20210513_S1/")),
  list(condition='glu8uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210513/20210513_curated/20210513_S2/", "./data_theo/20210708/20210708_curated/20210708_S0")),
  list(condition='glu16uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210513/20210513_curated/20210513_S3/")),
  list(condition='glu32uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210513/20210513_curated/20210513_S4/", "./data_theo/20210708/20210708_curated/20210708_S1")),
  list(condition='glu64uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210513/20210513_curated/20210513_S5/")),
  list(condition='glu128uM', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210708/20210708_curated/20210708_S2")),
  list(condition='glu002', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210708/20210708_curated/20210708_S4")),
  list(condition='glu02', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210708/20210708_curated/20210708_S6")),
  
  list(condition='lac002', duration=48*60, dt=180, medium='lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S1/", "./data_theo/20210305/20210305_curated/20210305_S1/")),
  list(condition='lac002_glu2uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S2/", "./data_theo/20210305/20210305_curated/20210305_S2/")),
  list(condition='lac002_glu4uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S3/", "./data_theo/20210305/20210305_curated/20210305_S3/")),
  list(condition='lac002_glu8uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S4/", "./data_theo/20210305/20210305_curated/20210305_S4/")),
  list(condition='lac002_glu16uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S5/", "./data_theo/20210305/20210305_curated/20210305_S5/", 
               "./data_theo/20210504/20210504_curated/20210504_S0/")),
  list(condition='lac002_glu32uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S6/", "./data_theo/20210305/20210305_curated/20210305_S6/",
               "./data_theo/20210504/20210504_curated/20210504_S1/")),
  list(condition='lac002_glu64uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210504/20210504_curated/20210504_S2/", "./data_theo/20210506/20210506_curated/20210506_S2/")),
  list(condition='lac002_glu128uM', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210504/20210504_curated/20210504_S3/", "./data_theo/20210506/20210506_curated/20210506_S3/", 
               "./data_theo/20210708/20210708_curated/20210708_S7")),
  list(condition='lac002_glu002', duration=48*60, dt=180, medium='glucose+lactose',
       paths=c("./data_theo/20210122/20210122_curated/20210122_S7/", "./data_theo/20210305/20210305_curated/20210305_S7/")),
  
  list(condition='glu2uM+IPTG', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210504/20210504_curated/20210504_S4/", "./data_theo/20210506/20210506_curated/20210506_S4/", 
               "./data_theo/20210513/20210513_curated/20210513_S6/")),
  list(condition='glu8uM+IPTG', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210504/20210504_curated/20210504_S5/", "./data_theo/20210506/20210506_curated/20210506_S5/")),
  list(condition='glu32uM+IPTG', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210504/20210504_curated/20210504_S6/", "./data_theo/20210506/20210506_curated/20210506_S6/")),
  list(condition='glu128uM+IPTG', duration=48*60, dt=180, medium='glucose',
       paths=c("./data_theo/20210506/20210506_curated/20210506_S7/", "./data_theo/20210513/20210513_curated/20210513_S7/"))
  
)

# LOAD MoMA DATA ####
# find raw data files from myconditions and store them in a dataframe
myfiles <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ which(names(.l) %in% c("condition", "paths"))] %>% 
           as.data.frame(stringsAsFactors=FALSE) ) %>% 
  do.call(rbind, .) %>% 
  rename(data_path=paths) %>%
  # for each path, find all files matched by the pattern .*\\d+\\.csv (e.g. *20151023.csv)
  group_by(condition, data_path) %>% 
  do((function(.df)
    # list.files(.df$data_path, ".*\\d+\\.csv", recursive=TRUE, full.names=TRUE) %>% 
    find.files(.df$data_path, "ExportedCellStats_*.csv") %>% 
      data.frame(path=., stringsAsFactors=FALSE) )(.))  

#  create condition_acq_times (describing acquisition times and temporal change of each condition) from myconditions
condition_acq_times <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ - which(names(.l) == "paths")] %>% as_tibble ) %>% 
  bind_rows() %>% 
  group_by(condition) %>% 
  mutate(m_start=cumsum(c(0, duration[-(length(duration))])) * 60,
         m_end=cumsum(duration) * 60 - 1e-5, duration=NULL,
         m_cycle=value_occurence_index(medium), 
  ) %>% 
  group_by(condition, m_start) %>% 
  mutate(data=list(data.frame(time_sec=seq(m_start, m_end, dt)) )) %>% unnest(cols=data) %>% 
  group_by(condition) %>% 
  mutate(frame=as.numeric(order(time_sec)-1))

# load perl scripts output to dataframes (using parallel dplyr)
myframes <- myfiles %>% 
  # process exported files on the cluster if required (otherwise return the list of paths)
  ungroup %>% 
  mutate(ppath=process_moma_data(path, .data2preproc=data2preproc, .frames_pl_script="get_size_and_fluo_multich.pl", #.skip=TRUE,
                                 .qsub_name="MMex_pl", .force=FALSE, .skip=TRUE) ) %>% 
  filter(!is.na(ppath)) %>% 
  # load perl scripts output to dataframes (in parallel, using multidplyr)
  group_by(condition, path) %>%
  partition(mycluster %>% cluster_copy(c("dl", "vertical_cutoff", "condition_acq_times"))) %>%
  do((function(.df){
    # browser()
    # print(.df$ppath)
    parse_frames_stats(.df$ppath)
  })(.)) %>%
  collect() %>%
  ungroup %>%
  extract(path, c("date", "pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  # propagate time and medium info
  left_join(condition_acq_times) %>%
  group_by(path) %>% 
  mutate(m_end=ifelse(m_end > max(time_sec), max(time_sec), m_end)) %>% 
  # fix end_type for pruned cells
  mutate(discard_start=(time_sec < 2*3600), length_um=length_pixel*dl) %>%
  group_by(path) %>%
  mutate(ndgt=compute_daughters_numbers(cid)) %>%
  mutate(end_type_moma=end_type,
         end_type=ifelse(ndgt==0, "lost", "weird"),
         end_type=ifelse(ndgt==2, "div", end_type)) %>%
  # remove frames after touching the exit
  group_by(id, .add=TRUE) %>%
  mutate(discard_top=which_touch_exit(vertical_top, vertical_cutoff)) %>%
  mutate(discard_top=ifelse(discard_start, FALSE, discard_top)) %>% # not in the preexpt step (2h)
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>% # update end_type to exit for cells which have touched the vertical cutoff
  # remove daughters of cells that touched the exit
  group_by(path) %>%
  mutate(discard_top=which_to_progeny(discard_top, cid)) %>%
  # append useful variables (per cell)
  mutate(
    cell=paste(date, pos, gl, id, sep='.'),
    ugen=paste(date, pos, gl, genealogy, sep='.'),
    gl_id=paste(date, pos, gl, sep='.'),
    vertical_center=(vertical_bottom + vertical_top)/2,
    mstep=paste(medium, m_cycle, sep='.'),
    strain='ASC662', strain=ifelse(condition=='mg1655', 'MG1655', strain),
    strain=ifelse(condition=='switch_long_lac_hiExpr', 'MG1655_pHi-GFP', strain),
    strain=ifelse(condition=='switch_∆lacA', 'AB460', strain),
  ) %>%
  group_by(date, pos, gl, id) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         b_rank=round(mean(total_cell_in_lane - cell_num_in_lane)),
         length_raw=(vertical_bottom-vertical_top)*dl,
         length_erik=length_um, length_um=length_raw)


# CONVERT FLUO UNITS ####
if (use_kaiser2018_params)
  autofluo_predict <- function(.h) .h * 422.8
myframes <- myframes %>% ungroup() %>% mutate(
  fluo_amplitude=fluo_amplitude * ifelse(date==20180712, 4, 1),
  fluo_amplitude=ifelse(date %in% c(20181008, 20181009) & fluo_amplitude==-1, NA, fluo_amplitude),
  fluo_amplitude=fluo_amplitude * ifelse(date==20181008 & pos %in% 0:4 & between(time_sec/3600, 8, 20-6/60), 5, 1),
  fluo_amplitude=fluo_amplitude * ifelse(date==20181009 & pos %in% c(0,2,4,6,8) & between(time_sec/3600, 8, 20-6/60), 5, 1),
  fluo_amplitude=fluo_amplitude * ifelse(date==20181024 & pos %in% c(0:2,4,6,8) & between(time_sec/3600, 8, 20-6/60), 5, 1),
  
  fluo_amplitude=ifelse(str_detect(condition, "(?:glu\\d+)|(?:lac\\d+)") & frame %% 3, NA, fluo_amplitude),
    
  fluogfp_amplitude = fluo_amplitude - autofluo_predict(length_um)
)


fp_per_oligomers <- 4 # lacZ is tetrameric
if (use_kaiser2018_params)
  fp_per_dn <- 0.0361 * fp_per_oligomers
myframes <- myframes %>% ungroup() %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  group_by(date, pos, gl, id)

# experiments to be discarded as identified by controls
discarded_dates <- c(
  # # slow growth in initial condition
  # 20151218, # switch_08h
  # 20161021, # switch_late
  # 20160526, # switch_12h_old (long gr tail + bad r2)
  # 20170108, # switch_late
  # 20170901, # switch_ramp 15'
  # 20180313, # switch_24h only 10GLs analysed
  # # other reasons
  # 20180123, # switch_lacIoe (weird lag distrib, but all longer than without overexpressing LacI)
  # 20180615,  # weird late switch control (all switch fast!)
  # 20180606  # weird late switch control (all switch fast!)
)


# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global envt (hence inheriting existing variables and keeping the newly created ones)...

rename_conds <- function (.str) {
# TODO: check NAs for all conditions
  .labels <- .str
  # .labels[.labels=='1'] <- 'switch:1'
  .labels <- str_replace(.labels, "switch_glcLac_lac", "glc+lac > lac")
  .labels <- str_replace(.labels, "switch_lactulose", "> lactulose")
  .labels <- str_replace(.labels, "switch_glycerol_", "> glyc ")
  # .labels <- str_replace(.labels, "glycerol", "glyc")
  .labels <- str_replace(.labels, "switch_0?(\\d+)h", "memory \\1h")
  .labels <- str_replace(.labels, "^switch_", "")
  .labels <- str_replace(.labels, "_", " ")
  return(.labels)
}
lac_lags_label <- expression(paste(italic('lac'), ' induction lag (min)'))


myscales <- list()
myplots <- list()
mytables <- list()

library(svglite)
knitr::opts_chunk$set(
  echo=FALSE, message=FALSE, warning=FALSE,
  dev="svglite"
)
# rmarkdown::clean_site()

# render control plots of each GC
# source('./src/MoM_lacDilution_GCplots.R')

# rmarkdown::render_site('./src/MoM_lacDilution_GFP_Estimation.Rmd')
# rmarkdown::render_site('./src/MoM_lacDilution_Constant_Envts.Rmd')
rmarkdown::render_site('./src/GCS_MoM_lac_Lags_Estimation.Rmd')

# DISCARD SOME DATASETS
rmarkdown::render_site('./src/GCS_MoM_lac_Controls.Rmd')

# rmarkdown::render_site('./src/index.Rmd')
# rmarkdown::render_site('./src/GCS_MoM_lac_Native.Rmd')
rmarkdown::render_site('./src/GCS_MoM_lac_Sensitivity.Rmd')


# RENDER ARTICLE FIGURES ####

knitr::opts_chunk$set(
  echo=FALSE, message=FALSE, warning=FALSE
  )

myfigs <- list()
source('./src/MoM_lacInduction_Figs.R')
source('./src/MoM_lacInduction_FigsSI.R')

# save.image(".RData")
