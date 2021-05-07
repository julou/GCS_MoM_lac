# # SAVE ENVIRONMENT (but `pls`)
# myvars <- ls(all.names = TRUE)
# save(list=myvars[which(myvars!='pls')], file=".RData", envir=.GlobalEnv)

# myframes %>%
#   filter((date=='20150703' & pos==0 & gl==9)) %>%
#   filter(!discard_top, vertical_top>vertical_cutoff) %>%
#   mutate(time_sec=time_sec-2*3600) %>%
#   plot_faceted_var_tracks(.var_col='length_um', .show_cellid=TRUE, .log=TRUE) +
#   # mask early frames (requires a dummy df!)
#   geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
#   scale_x_hours(4) +
#   scale_y_continuous(breaks=2:4, trans='log2') +
#   labs(y='length (µm)')

# Plot overall experiment
myframes %>% 
  # filter(date==20170901) %>% 
  # create a list of filenames for all plots
  group_by(condition, strain, date) %>% 
  slice(1) %>% 
  select(condition, strain, date) %>% 
  na.omit %>% 
  mutate(., basename=paste(condition, date, sep="_")) %>% 
  left_join(., expand.grid(condition=unique(.[['condition']]), var=c("length_um", "gfp_nb", "gfp_conc"), stringsAsFactors=FALSE)) %>% 
  left_join(data.frame(var=c("length_um", "gfp_nb", "gfp_conc"), log=c(TRUE, TRUE, FALSE))) %>% 
  mutate(filename=file.path(".", "plots", "GCplots", sprintf("%s__%s.pdf", basename, var))) %>% 
  # for each file, create it if missing (partition() doesn't work yet)
  group_by_(., .dots=names(.)) %>% # grouping by all variables
  # partition_(., groups=names(.), cluster=mycluster %>% cluster_copy(c("myframes", "vertical_cutoff")) %>% 
  #              cluster_copy("plot_faceted_var_tracks")) %>% # grouping by all variables
  do((function(.dff, .skip=TRUE){
    if (file.exists(.dff$filename) && .skip) return(data.frame())
    # browser()
    .pls <- myframes %>% 
      ungroup %>% semi_join(as.data.frame(.dff)) %>% 
      mutate(gfp_conc=gfp_nb/length_um) %>% 
      # top_n(1e4) %>%
      group_by(condition, strain, date, pos, gl) %>%
      do(pll=(function(.dfgl, .b_rank_max=6){
        # browser()
        custom_labels <- function (.str) {
          .labels <- paste('rank:', .str)
          .labels[.labels=='rank: -1'] <- 'all'
          .labels[.labels==paste0('rank: ', .b_rank_max)] <- paste0('rank: >= ', .b_rank_max)
          return(.labels)
        }
        
        # if (nrow(filter(.dfgl, !discard_top, vertical_top>vertical_cutoff)) == 0) return(list())
        .dfgl <- .dfgl %>% ungroup %>% filter(!discard_top, vertical_top>vertical_cutoff) %>% 
          mutate(time_sec=time_sec-2*3600, b_rank=ifelse(b_rank>.b_rank_max, .b_rank_max, b_rank))
        if (nrow(.dfgl) == 0) return(list())
        if (.dff$var %in% c("gfp_nb", "gfp_conc")) .dfgl <- filter(.dfgl, !is.na(gfp_nb))
          
        .pl <- .dfgl %>% 
          plot_faceted_var_tracks(.var_col=.dff$var, .show_all=TRUE, .show_cellid=TRUE, .log=.dff$log, .facet_labeller=custom_labels, .dt=first(.dfgl$dt)) +
          # mask early frames (requires a dummy df!)
          geom_rect(xmin=-Inf, xmax=0, ymin=ifelse(.dff$log, 0, -Inf), ymax=Inf, fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
          scale_x_hours(4) 
        # labs(title=sprintf("%s(%s) %s.%02d.%02d", .dff$strain, .dff$condition, .dff$date, first(.dfgl$pos), first(.dfgl$gl) ))
        
        if (.dff$var=='length_um') 
          .pl <- .pl + scale_y_continuous(breaks=2:4, trans='log2') + labs(y='length (µm)')
        
        return(.pl)
      })(.))
    
    # Print plots to pdf
    pdf(.dff$filename, width=12, height=10)
    for (.i in 1:nrow(.pls)) {
      try(plot(.pls[[.i, 'pll']] + 
                 labs(title=sprintf("%s(%s) %s.%02d.%02d", .pls[[.i, "strain"]], .pls[[.i, "condition"]], 
                                    .pls[[.i, "date"]], .pls[[.i, "pos"]], .pls[[.i, "gl"]]))) )
    }
    dev.off()
    
    return(data.frame())
  })(.))


