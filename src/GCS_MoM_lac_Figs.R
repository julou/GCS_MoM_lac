myfigs <- list() 

### ### ### ###
#### MODELING PLOTS ####

(myplots[['lac_model_induction']] <-
   tibble(path=list.files(here("material", "Erik_models"), "lac_induc.*", full.names = TRUE)) %>% 
   mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
   unnest(data) %>% 
   extract(path, c('curve', 'branch'), ".*/lac_induc(\\d)_(.+)\\.txt") %>% 
   ggplot() +
   geom_line(aes(exp(X1), exp(X2), lty=branch, col=curve)) +
   scale_x_log10(breaks=c(1e-4, 1e-3, 1e-2, .1), expand=c(0, 0),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
   scale_y_log10(breaks=c(10, 1e3, 1e5),
                 # breaks = trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
   expand_limits(y=c(2.5, 1e5)) +
   labs(x=expression(paste('inducer level ', italic('b'))), y='LacY molecules') +
   scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dotted')) +
   scale_color_manual(values=c('1'='black', '2'='grey60')) +
   theme(legend.position = 'none') +
   NULL
)

(myplots[['lac_model_const_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "lac_phasediagram_cons.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/lac_phasediagram_cons_(.+)\\.txt") %>% 
    ggplot() +
    geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                 data=~bind_rows(., 
                                 tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                 tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')),
                 )) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_vline(xintercept = log(2)/0.5, lty='dashed') +
    geom_vline(xintercept = log(2)/0.5*10, lty='dashed', col='grey60') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-6, 1e-3, 1), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(y=c(5e-7, 2)) +
    labs(x='doubling time (h)', y=expression(paste('inducer level ', italic('b')))) +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['lac_model_crp_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "lac_phasediagram_w_crp.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/lac_phasediagram_w_crp_(.+)\\.txt") %>% 
    ggplot() +
    geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                 data=~bind_rows(., tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_vline(xintercept = log(2)/0.5, lty='dashed') +
    geom_vline(xintercept = log(2)/0.5*10, lty='dashed', col='grey60') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-6, 1e-3, 1), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.3, y=5e-7) +
    labs(x='doubling time (h)', y=expression(paste('inducer level ', italic('b')))) +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['autoactiv_model_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "autoactivating_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/autoactivating_phasediagram_(.+)\\.txt") %>% 
    ggplot() +
    geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                 data=~bind_rows(., tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                 tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) )) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_hline(yintercept = 0.061, lty='dashed') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    # annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-3, 1e-1, 10), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.35, y=c(2e-4, 20)) +
    labs(x='doubling time (h)', y='basal production rate       ') +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['autoactiv_model_induction']] <-
    tibble(path=list.files(here("material", "Erik_models"), "autoactivating_induc_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('branch'), ".*/autoactivating_induc_(.+)\\.txt") %>% 
    ggplot() +
    geom_line(aes(exp(X1), exp(X2), lty=branch)) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-1, 1e1, 1e3),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x='doubling time (h)', y=expression(paste('TF concentration ', italic('x'), '            '))) +
    scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dotted')) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['signal_dep_decay_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "signal_coupled_decay_pd_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest() %>% 
    extract(path, c('curve'), ".*/signal_coupled_decay_pd_(.+)\\.txt") %>% 
    # with(range(X1))
    ggplot() +
    geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                 data=~bind_rows(., tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                 tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) )) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_hline(yintercept = 100, lty='dashed', col='gray50') +
    geom_hline(yintercept = 1e4, lty='dashed') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    # annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    # scale_x_log10(limits=c(.5, 1e3), breaks=c(1, 24, 500), expand=c(0, 0)) +
    scale_x_log10(limits=c(.5, 1000), breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(limits=c(.4, 1e5), breaks=c(1, 1e2, 1e4), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    # labs(x='doubling time (h)', y=expression(paste('signal strength ', italic('s'), '/', italic('s[0]')))) +
    labs(x='doubling time (h)', y=expression(paste('signal strength ', s / s[0], '     ' ))) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['signal_dep_decay_induction']] <-
    tibble(path=list.files(here("material", "Erik_models"), "competence_induc_s100.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    extract(path, c('signal', 'branch'), ".*/competence_induc_s(\\d+)_(.+)\\.txt") %>% 
    unnest(data) %>% 
    ggplot() +
    geom_line(aes(exp(X1), exp(X2), lty=branch, col=signal, group=interaction(signal, branch))) +
    # scale_x_log10(limits=c(.5, 1e3), breaks=c(1, 24, 480), expand=c(0, 0)) +
    scale_x_log10(limits=c(.5, 1000), breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-1, 1e1, 1e3),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    labs(x='doubling time (h)', y=expression(paste('TF concentration ', italic('x'), '            '))) +
    scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dotted')) +
    scale_color_manual(values=c('100'='grey50', '10000'='black')) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['twocmps_model_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "twocomp_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/twocomp_phasediagram_(.+)\\.txt") %>% 
    # with(range(X1))
    ggplot() +
    geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                 data=~bind_rows(., tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                 tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) )) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_vline(xintercept = 2, lty='dashed', col='gray50') +
    geom_vline(xintercept = 20, lty='dashed') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    # annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    # scale_x_log10(limits=c(.5, 500), breaks=c(1, 24, 480), expand=c(0, 0)) +
    scale_x_log10(limits=c(.5, 500), breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-6, .01, 1e2), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.5, y=c(2e-8, 7)) +
    # labs(x='doubling time (h)', y=expression(paste('signal strength ', italic('s'), '/', italic('s[0]')))) +
    labs(x='doubling time (h)', y=expression(paste('signal strength ', s / s[0], '     ' ))) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['twocmps_model_induction']] <-
    tibble(path=list.files(here("material", "Erik_models"), "twocomp_induc_td2.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    extract(path, c('dt', 'branch'), ".*/twocomp_induc_td(\\d+)_(.+)\\.txt") %>% 
    unnest(data) %>% 
    ggplot() +
    geom_line(aes(exp(X1), exp(X2), lty=branch, col=dt, group=interaction(dt, branch))) +
    scale_x_log10(breaks=c(.01, .1, 1, 10), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-1, 1e1, 1e3),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(y=c(1e-1, 1e3)) +
    labs(x=expression(paste('signal strength ', s / s[0], '     ' )), y=expression(paste('TF concentration ', italic('x'), '            '))) +
    scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dotted')) +
    scale_color_manual(values=c('2'='grey50', '20'='black')) +
    theme(legend.position = 'none') +
    NULL
)

(myplots[['critic_conc_monod']] <- (function() {
  lamstar <- 1.2
  mux <- 1/48
  muy <- 1/48
  ggplot() +
    # stat_function(fun=~ (.+mux)/(lamstar-.), position='jitter') +
    stat_function(aes(ymin=after_stat(log10(y)), ymax=Inf, fill="lower"), n=501,
                  fun=~ 7.9^(-1/2) * (.+mux)/(lamstar-.), geom='ribbon') +
    stat_function(aes(y=after_stat(log10(y)), col="lower"), n=501,
                  fun=~ 7.9^(-1/2) * (.+mux)/(lamstar-.)) +
    stat_function(aes(ymin=after_stat(log10(y)), ymax=Inf, fill="upper"), n=501,
                  fun=~ 7.9^(1/2) * (.+mux)/(lamstar-.), geom='ribbon') +
    stat_function(aes(y=after_stat(log10(y)), col="upper"), n=501,
                  fun=~ 7.9^(1/2) * (.+mux)/(lamstar-.)) +
    stat_function(aes(y=after_stat(log10(y))), fun=~ ./(lamstar-.), lty='dotted') +
    annotate("text", x=1.2, y=-2, label='uninduced', hjust=1.5, vjust=0) +
    annotate("text", x=0, y=Inf, label='induced', hjust=-.1, vjust=2) +
    coord_cartesian(xlim=c(-.0, 1.2), ylim=c(-2.6, 3), expand = FALSE) +
    # coord_cartesian(xlim=c(-0.01, 1.2), expand = FALSE) +
    scale_x_continuous(limits=c(-.01, 1.2), breaks=scales::breaks_pretty(n=4)) +
    scale_y_continuous(breaks=seq(-2, 2, 2), labels = 10^seq(-2, 2, 2)) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    labs(x="growth rate (dbl/h)", y="critical nutrient concentration c/c0") +
    guides(col='none', fill='none') +
    NULL
})() )


# ### ### ### ###
# #### MILLER ASSAY ####
(function() {
  # load TMG plots into global env
  load(list.files("data", "GCS_Miller_lac_plots_\\w+\\.RData", full.names=TRUE) ) # do NOT use load() after %>%
  # load('data/TMG_sensitivity_Miller_plots.RData')
  myplots_miller <<- myplots
  # TODO: load from RData file to global env
  mysugars <<- c('0.2% arabinose+CA'='M9+0.2ara+CA', '0.2% arabinose'='M9+0.2ara', '0.2% glycerol'='M9+0.2gly', 
                '0.2% mannose'='M9+0.2man', '0.2% pyruvate'='M9+0.2pyr',
                '0.2% ribose'='M9+0.2rib', '0.2% succinate'='M9+0.2suc', '0.4% rhamnose'='M9+0.4rhm')
  scale_color_sugars <<- function(...) 
    ggplot2::scale_color_manual(..., values = c(qual_cols, "gray30"), na.value='gray70', 
                                breaks = mysugars, labels = names(mysugars), 
                                limits = function(.l) intersect(mysugars, .l) )
  
})()


### ### ### ###
#### FIG Models ####
(myfigs[['GCS_models']] <- plot_grid(
    NULL, 
    myplots[['autoactiv_model_phdiag']], 
    myplots[['autoactiv_model_induction']], 
    NULL, 
    myplots[['signal_dep_decay_phdiag']], 
    myplots[['signal_dep_decay_induction']], 
    NULL,
    myplots[['twocmps_model_phdiag']], 
    myplots[['twocmps_model_induction']], 
    nrow=3, byrow=FALSE, align='vh', labels=LETTERS[c(1, 4, 7, 2, 5, 8, 3, 6, 9)], rel_heights=c(1.2, 1.2, 1)
) )

save_plot(here("plots", "figs", "GCS_MoM_lac_fig_models.pdf"), myfigs[['GCS_models']],
          base_height=NULL, base_width=4.75 * 14/7, # 2 cols
          base_asp = 1.5
)


### ### ### ###
#### FIG GCS lac operon ####
(myfigs[['GCS_lac']] <- plot_grid(
  # LEFT COL
  # plot_grid(
  #   NULL,
  #   myplots[['lac_model_const_phdiag']],
  #   # ggdraw( myplots[['lac_model_const_phdiag']] ) +
  #   #   draw_plot( myplots[['lac_model_induction']] , .45, .45, .5, .5),
  #   ncol=1, rel_heights = c(1, 1.1), labels=c("A", "B")
  # ),
  plot_grid(
    NULL,
    myplots[['lac_model_const_phdiag']] +
      labs(y=expression(paste('inducer level ', italic('b'), '      ')) ),
    myplots[['lac_model_induction']] +
      labs(y="LacY molecules         "),
    ncol=1, rel_heights = c(1, 1, 1), labels=c("AUTO")
  ),
  
  # MID COL
  plot_grid(
    NULL,
    
    plot_grid(
      myplots[['lacl_boxpl_subset']](igr, 
                                     condition %in% c("switch_lactulose_TMG20", "switch_glycerol_TMG20"),
                                     time_bin %in% c("[-1,0]", "(0,1]", "(7,8]"),
                                     .simplify_time_labels=TRUE) +
        geom_rect(xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf, fill='black', alpha=.03) +
        scale_colour_manual(values=c("non induced"="gray25", "induced"=qual_cols[3], "non growing"=qual_cols[1], "growing"=qual_cols[2]),
                            breaks=c("non induced", "induced", "growing", "non growing"),
                            guide=guide_legend(ncol=2)) +
        scale_y_continuous(breaks=c(0, 0.5, 1)) +
        coord_cartesian(ylim=c(-.1, 1.1)) +
        labs(x=NULL, y="growth  \nrate (dbl/h)  ") +
        theme(strip.text.x = element_blank(), legend.position = "none", 
              axis.text.x = element_blank(), # axis.title.x = element_blank(), 
              panel.spacing.x = unit(0.5, "lines"),
              # axis.title.y = element_text(margin=margin(t=30)), 
              axis.text.y = element_text(margin = margin(r=3, l=-8)),
              plot.margin = margin(7, 7, 2, 7),
        ) +
        NULL,
      
      myplots[['lacl_boxpl_subset']](gfp_c,
                                     condition %in% c("switch_lactulose_TMG20", "switch_glycerol_TMG20"),
                                     time_bin %in% c("[-1,0]", "(0,1]", "(7,8]"),
                                     .simplify_time_labels=TRUE) +
        geom_rect(xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf, fill='black', alpha=.03) +
        scale_colour_manual(values=c("non induced"="gray25", "induced"=qual_cols[3], "non growing"=qual_cols[1], "growing"=qual_cols[2]),
                            limits=c("non induced", "induced", "growing", "non growing"), 
                            guide=guide_legend(ncol=2)
        ) +
        coord_cartesian(ylim=c(NA, 1400)) +
        labs(y="[LacZ-GFP]          \n(FP/µm)           ", col="") +
        theme_cowplot_legend_inset() +
        theme(strip.text.x = element_blank(), legend.position = "bottom", 
              panel.spacing.x = unit(0.5, "lines"),
              # axis.title.y = element_text(margin=margin(t=50)),
              axis.text.y = element_text(margin = margin(r=3, l=-8)),
              
              legend.margin = margin(t=-12, b=-4),
              plot.margin = margin(0, 7, 5, 7),
        ) +
        NULL,
      
      ncol=1, rel_heights = c(.4, .6), align = "v" ), 
    
    ncol=1, rel_heights = c(1, 1.1), labels = c("D", "")
  ),
  
  # RIGHT COL
  plot_grid(
    myplots_miller[['GCS_sugars_indct_subset']] +
      # coord_cartesian(xlim=c(.2, NA), ylim=c(.5, 1000)) +
      theme(legend.position = 'none') +
      labs(title=" ", y=expression('p'['lac']*' activity (MU/h)    ')) +
      NULL,
    myplots_miller[['GCS_sugars']]() + 
      guides(col=guide_legend(title=NULL, direction = "vertical", ncol=2, )) +
      theme_cowplot_legend_inset() +
      theme(legend.position = 'top', #legend.box.spacing=unit(-1, "line"),
            legend.box.margin=margin(t=-8, l=-48), # shift the whole legend horizontally
            # legend.margin=margin(t=-20),
      ) +
      labs(y='critical [TMG] (µM)  ') +
      NULL,
    ncol=1, rel_heights = c(1, 1.45), labels = c("E", "F"), align="v", axis="l"
  ),
  
  nrow = 1, rel_widths = c(.8, 1.05, 1)
) )

save_plot(here("plots", "figs", "GCS_MoM_lac_fig_lac.pdf"), myfigs[['GCS_lac']],
          base_height=NULL, base_width=4.75 * 14/7, # 2 cols
          base_asp = 1.8 )


### ### ### ###
#### FIG native + cm ####
(myplots[['Cline_You2013']] <- 
   bind_rows(
     tribble(
       # Table S1 from You, et al. 2013
       ~c_source, ~gr, ~gr_se, ~LacZ, ~LacZ_se,
       "60 mM acetate", 0.37, 0.00, 20.8, 0.3,
       "20 mM arabinose", 0.40, 0.01, 22.8, 0.1,
       "20 mM mannose", 0.41, 0.00, 20.4, 1.8,
       "15 mM succinate", 0.46, 0.01, 20.4, 0.5,
       "20 mM sorbitol", 0.46, 0.01, 18.4, 0.2,
       "20 mM pyruvate", 0.61, 0.02, 15.6, 1.7,
       "20 mM fructose", 0.61, 0.02, 17.3, 0.1,
       "0.4% (v/v) glycerol", 0.63, 0.01, 16.4, 0.2,
       "0.2% (w/v) maltose", 0.67, 0.00, 15.8, 0.3,
       "0.4% (w/v) glucose", 0.85, 0.00, 8.23, 0.21,
       "20 mM gluconate", 0.88, 0.04, 9.45, 0.99,
       "0.2% (w/v) lactose", 0.98, 0.01, 5.80, 0.09,
       "10 mM glucose-6P+10 mM gluconate", 1.09, 0.01, 1.87, 0.14,
     ) %>% mutate(type="C source"), 
     tribble(
       # Table S2 from You, et al. 2013
       ~mba , ~gr, ~gr_se, ~LacZ, ~LacZ_se,
       0, 0.39, 0.01, 18.8, 0.9, 
       12.5, 0.48, 0.01, 17.2, 0.8, 
       25, 0.57, 0.01, 14.9, 0.7, 
       50, 0.63, 0.01, 13.2, 0.6, 
       100, 0.70, 0.01, 11.4, 0.5, 
       200, 0.78, 0.02, 12.0, 0.6, 
       500, 0.84, 0.02, 9.58, 0.44
     ) %>% mutate(type='titr. LacY')
   ) %>% 
   mutate(LacZ=1e3*LacZ, LacZ_se=1e3*LacZ_se) %>% 
   # filter()
   ggplot(aes(gr, LacZ, col=type)) +
   # stat_smooth(aes(col=type), method='lm', se=FALSE, fullrange=TRUE) +
   geom_errorbar(aes(ymin=LacZ-LacZ_se, ymax=LacZ+LacZ_se, width=0), data=~filter(., type=='C source')) +
   geom_errorbarh(aes(xmin=gr-gr_se, xmax=gr+gr_se, height=0), data=~filter(., type=='C source')) +
   geom_point(size=2, data=~filter(., type=='C source')) +
   stat_smooth(method='lm', se=FALSE, fullrange=TRUE, col='gray20') +
   geom_errorbar(aes(ymin=LacZ-LacZ_se, ymax=LacZ+LacZ_se, width=0), data=~filter(., type=='titr. LacY')) +
   geom_errorbarh(aes(xmin=gr-gr_se, xmax=gr+gr_se, height=0), data=~filter(., type=='titr. LacY')) +
   geom_point(size=2, data=~filter(., type=='titr. LacY')) +
   scale_colour_grey(start=0, end=.8, limits=c('titr. LacY', 'C source'), guide='none') +
   expand_limits(x=c(0, 1.2), y=33) +
   labs(x='growth rate (dbl/h)', y='max [LacZ] (MU)') +
   NULL)


# ggplot() +
#   stat_function(aes(col='CRP-like'), fun=function(.x) 1-.x) +
#   stat_function(aes(col='other'), fun=function(.x) exp(-.x*5)) +
#   xlim(0, 1) +
#   labs(x='growth rate', y='max expression', col='regulation') +
#   theme(axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.position = c(1,1), legend.justification = c(1,1)) +
#   NULL

(myfigs[[2]] <- plot_grid(
  NULL, 
  NULL,
  
  myplots[['Cline_You2013']] +
    scale_x_continuous(breaks=c(0, 0.5, 1)) +
    coord_cartesian(xlim=c(0, 1.25), ylim=c(0, NA)) +
    expand_limits(x=c(-.1, 1.4)) +
    labs(y='max [LacZ] (MU)      ') +
    NULL,
  myplots[['lac_model_crp_phdiag']],
  
  NULL,
  NULL,

  myplots_miller[['GCS_cm_maxindct']](strain=='MG1655') +
    scale_x_continuous(breaks=c(0, 0.5, 1)) +
    coord_cartesian(xlim=c(0, 1.25), ylim=c(0, NA)) +
    expand_limits(x=c(-.1, 1.4)) +
    labs(y='max [LacZ] (MU)      ') +
    theme(legend.position = 'none') +
    NULL,
  myplots_miller[['GCS_cm']](strain=='MG1655') +
    guides(shape='none') +
    scale_color_viridis_c(breaks=c(0,4,8)) +
    theme_cowplot_legend_inset() +
    labs(y="critical [TMG] (µM)      ", col="[cam]\n(µM)") +
    theme(legend.position = "right", legend.key.height = unit(10, "pt")) +
    guides(col=guide_colourbar(direction = "horizontal", title.position = "left", barwidth=unit(58, 'pt'))) +
    theme(legend.position = c(0.03,0.03), legend.justification = c(0,0)) +
    NULL,
  
  ncol=2, labels=c('', '', 'A', 'B', '', '', 'C', 'D'),
  rel_widths=c(1, 1.1), rel_heights=c(.1, 1, .1, 1), 
  align='hv'
) )

save_plot(here("plots", "figs", "GCS_MoM_lac_fig_cm.pdf"), myfigs[[2]],
          base_height=NULL, base_width=2.25 * 14/7, # 1 col
          base_asp = 1
)


### ### ### ###
#### FIG GCS with regulation ####
(myfigs[['GCS_regul']] <- plot_grid(
  # plots row
  plot_grid(
    myplots[['critic_conc_monod']] +
      # labs(y=expression(paste("critical nutrient concentration ", c/c[0])) ) +
      labs(y="critical nutrient\nconcentration c/c0") +
      # labs(y="critical nutrient<br/>concentration c/c<sub>0</sub>") +
      # theme(axis.title.x = element_markdown()) +
      NULL,
    myplots[['sugarmix_crp']] +
      expand_limits(x=c(-.1, 1.3)) +
      coord_cartesian(xlim=c(0, 1.15), ylim=c(0, NA)) +
      # theme(legend.position = 'right') +
      theme(legend.position = c(1, 1), legend.justification = c(1, 1)) +
      NULL,
    ncol=1, rel_heights = c(1, 1), labels = c("A", "B")),
  
  myplots[['sugarmix_induction']](.xbreaks = 2 * 10^(0:4)),
  
  plot_grid(
    myplots_miller[['GCS_cm_maxindct']](strain=='MG1655') +
      scale_x_continuous(breaks=c(0, 0.5, 1)) +
      coord_cartesian(xlim=c(0, 1.25), ylim=c(0, NA)) +
      expand_limits(x=c(-.1, 1.4)) +
      labs(y='max [LacZ] (MU)      ') +
      theme(legend.position = 'none') +
      NULL,
    myplots_miller[['GCS_cm']](strain=='MG1655', .with_sugar_data=TRUE) +
      # geom_point(aes(x2, y2), col="gray85", data=
      #              myplots_miller[['GCS_sugars']]() %>% layer_data(2) %>% mutate(x2=exp(x), y2=exp(y)) ) +
      guides(shape='none') +
      scale_color_viridis_c(breaks=c(0,4,8)) +
      theme_cowplot_legend_inset() +
      labs(y="critical [TMG] (µM)      ", col="[cam]\n(µM)") +
      theme(legend.position = "right", legend.key.height = unit(10, "pt")) +
      guides(col=guide_colourbar(direction = "horizontal", title.position = "left", barwidth=unit(58, 'pt'))) +
      theme(legend.position = 'top') +
      # theme(legend.position = c(0.03,0.03), legend.justification = c(0,0)) +
      NULL,
    ncol=1, rel_heights = c(1, 1.2), labels = c("D", "E")), align='v',
  
  nrow=1, rel_widths = c(1.1, 1.5, 1), labels=c("", "C", "")
) )
save_plot(here("plots", "figs", "GCS_MoM_lac_fig_regul.pdf"), myfigs[['GCS_regul']],
          base_height=NULL, base_width=4.75 * 14/7, # 2 cols
          base_asp = 1.8
)


