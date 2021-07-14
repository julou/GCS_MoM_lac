
### ### ### ###
#### FIG 1 ####

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
   labs(x='doubling time (h)', y=expression(paste('operon expression ', italic('x')))) +
   scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dotted')) +
   theme(legend.position = 'none') +
   NULL
)

(myplots[['autoactiv_model_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "autoactivating_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/autoactivating_phasediagram_(.+)\\.txt") %>% 
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                   tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) ))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_hline(yintercept = 0.014269, lty='dashed') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    # annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-4, 1e-2, 1), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.35, y=c(7e-5, 7)) +
    labs(x='doubling time (h)', y='promoter activity') +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)


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

(myplots[['lac_model_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "lac_phase.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/lac_phasediagram_(.+)\\.txt") %>% 
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper'))))) +
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

(myplots[['twocmps_model_phdiag']] <-
    tibble(path=list.files(here("material", "Erik_models"), "twocomp_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest(data) %>% 
    extract(path, c('curve'), ".*/twocomp_phasediagram_(.+)\\.txt") %>% 
    # with(range(X1))
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                   tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) ))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    # geom_hline(yintercept = 0.014269, lty='dashed') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    # annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 24, 1000), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-3, 1, 1e3), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.15, y=c(2e-6, 7)) +
    # labs(x='doubling time (h)', y=expression(paste('signal strength ', italic('s'), '/', italic('s[0]')))) +
    labs(x='doubling time (h)', y=expression(paste('signal strength ', s / s[0] ))) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

(myfigs[[1]] <- plot_grid(
  # plots row
  plot_grid(
    NULL, NULL, NULL,
    nrow=1, labels=c("A", "C", "E"), rel_widths=c(1, 1, 1)
  ),
  plot_grid(
    myplots[['autoactiv_model_phdiag']], 
    myplots[['lac_model_phdiag']], 
    myplots[['twocmps_model_phdiag']], 
    nrow=1, labels=c("B", "D", "F"), rel_widths=c(1, 1, 1)
  ),
  ncol=1, rel_heights=c(1, 1)
) )

save_plot(here("plots", "figs", "GCS_MoM_lac_fig1.pdf"), myfigs[[1]],
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_asp = 1.95
)

### ### ### ###
#### FIG 2 ####
(function() {
  # load TMG plots into global env
  load('data/TMG_sensitivity_Miller_plots.RData')
  myplots_tmg <<- myplots
})()

(myfigs[[2]] <- plot_grid(
  # LEFT COL
  plot_grid(
    myplots_tmg[['GCS_sugars_indct_subset']] +
      theme(legend.position = 'none') +
      # theme(legend.position = c(1, 0), legend.justification = c(1, 0)) +
      # theme_cowplot_legend_inset() +
      # coord_cartesian(xlim=c(.2, NA), ylim=c(.5, 1000)) +
      labs(title=" ", y=expression('p'['lac']*' activity (MU/h)        ')) +
      NULL,
    myplots_tmg[['GCS_sugars_subset']] + 
      # scale_color_sugars(guide=guide_legend(reverse=TRUE)) +
      theme_cowplot_legend_inset() +
      theme(legend.position = 'right') +
      labs(y='critical [TMG] (µM)    ') +
      NULL,
    ncol=1, rel_heights = c(0.8, 1), labels = c("A", "B"), align="v", axis="l"
  ),
  
  # RIGHT COL
  plot_grid(
    NULL,
    
    plot_grid(
      myplots[['lacl_gly_mainfig']](igr, 
                                    condition %in% c("switch_lactulose_TMG20", "switch_glycerol_TMG20"),
                                    time_bin %in% c("[-1,0]", "(0,1]", "(7,8]")) +
        geom_rect(xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf, fill='black', alpha=.03) +
        scale_colour_manual(values=c("non induced"="gray25", "induced"=qual_cols[3], "non growing"=qual_cols[1], "growing"=qual_cols[2]),
                            breaks=c("non induced", "induced", "growing", "non growing"),
                            guide=guide_legend(ncol=2)) +
        scale_y_continuous(breaks=c(0, 0.5, 1)) +
        coord_cartesian(ylim=c(-.1, 1.1)) +
        labs(x=NULL, y="growth  \nrate (dbl/h)  ") +
        theme(strip.text.x = element_blank(), legend.position = "none", 
              axis.text.x = element_blank(), # axis.title.x = element_blank(), 
              panel.spacing.x = unit(2, "lines"),
              # axis.title.y = element_text(margin=margin(t=30)), 
              plot.margin = margin(7, 7, 2, 7),
        ) +
        NULL,
      
      myplots[['lacl_gly_mainfig']](gfp_c,
                                    condition %in% c("switch_lactulose_TMG20", "switch_glycerol_TMG20"),
                                    time_bin %in% c("[-1,0]", "(0,1]", "(7,8]")) +
        geom_rect(xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf, fill='black', alpha=.03) +
        scale_colour_manual(values=c("non induced"="gray25", "induced"=qual_cols[3], "non growing"=qual_cols[1], "growing"=qual_cols[2]),
                            limits=c("non induced", "induced", "growing", "non growing"), 
                            # guide=guide_legend(ncol=2)
                            ) +
        coord_cartesian(ylim=c(NA, 1400)) +
        labs(y="[LacZ-GFP]          \n(/µm)          ") +
        theme_cowplot_legend_inset() +
        theme(strip.text.x = element_blank(), legend.position = "bottom", 
              panel.spacing.x = unit(2, "lines"),
              # axis.title.y = element_text(margin=margin(t=50)),
              
              legend.margin = margin(t=-12, b=-4, r=50),
              plot.margin = margin(0, 7, 5, 7),
        ) +
        NULL,
      
      ncol=1, rel_heights = c(.4, .6), align = "v" ), 
    
    ncol=1, rel_heights = c(1, 1), labels = c("C", "D")
  ),
  nrow = 1, rel_widths = c(1, 1.1)
) )

save_plot(here("plots", "figs", "GCS_MoM_lac_fig2.pdf"), myfigs[[2]],
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_asp = 1.95 )


#####
# myfigs[[3]] <- 
#   plot_grid(
#     myplots[['naive_arrest_frac']],
#     myplots[['lacl_arrest_induction_fracs']],
#     NULL,
#     labels="AUTO", ncol=1, rel_heights=c(0.45, 1, 1)
#   )
# 
# save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]],
#           base_height=NULL, base_width=2.25 * 14/8, # 1 col
#           base_aspect_ratio = 1/2.2
# )

############
(myfigs[[3]] <- function() { # local envt
  pdftools_installed <- require(pdftools)
  plot_grid(
    plot_grid(
      plot_grid(
        plot_grid(
          myplots[['lac_model_induction']] + theme(axis.title.y = element_text(hjust=0) ),
          myplots[['lac_model_phase']] + theme(axis.title.y = element_text(hjust=0) ),
          nrow=1, labels=c('A', 'B'), rel_widths=c(1, 1), align='h'),
        myplots[['TMG-Miller-induction']] +
          theme_cowplot_legend_inset() +
          theme(axis.title.y = element_text(hjust=0),
                #plot.margin = margin(t=10, b=10, l=6, r=10),
                legend.title = element_blank(),
          ) +
          NULL,
        ncol=1, labels=c('', 'C'), rel_heights=c(1, 1), align='v' ),
      myplots[['TMG-Miller-conc-dt']] +
        theme(#plot.margin = margin(t=10, b=10, l=6, r=10)
        ) +
        NULL,
      nrow=1, labels=c('', 'D'), rel_widths=c(2, 1)),
    # caption row
    plot_grid(NULL, get_legend(myplots[['lacl_gr_hist']]() + guides(col='legend') +
                                 scale_colour_discrete(name='nutrient', breaks=c('switch_glycerol_TMG20', 'switch_lactulose_TMG20'),
                                                       labels=c('glycerol', 'lactulose')) +
                                 theme_cowplot_legend_inset() +
                                 theme(legend.position = 'top', legend.justification = c(0, 0), 
                                       legend.box = 'vertical', legend.box.just = 'left', 
                                       legend.spacing = unit(0, 'mm'), legend.box.spacing = unit(0, 'mm'),
                                       legend.margin=margin(), legend.box.margin = margin(t=50),
                                 )),
              nrow=1, rel_widths=c(1, 1.9)), 
    # lactulose row
    plot_grid(
      if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "MoM_lacDilution_fig3_cartoon.ai.pdf"), pages=1), scale=1.0) + 
        theme(plot.margin = margin(t=4)),
      if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "MoM_lacDilution_fig3_cartoon.ai.pdf"), pages=2), scale=1.0) + 
        theme(plot.margin = margin(t=4)),
      # if (!pdftools_installed) NULL else ggdraw() + draw_image(here("material", "montage_TMG_glyc_lacl.jpg"), y=-0.0, scale=0.8),
      myplots[['lacl_gr_hist']]() +
        theme(plot.margin = margin(t=30),
              legend.position = 'none',
              strip.background = element_blank(),
              strip.text.y = element_blank()) +
        NULL, 
      myplots[['lacl_lacGFP_hist']]() +
        theme(plot.margin = margin(t=30, l=10, r=10),
              legend.position = 'none',
              axis.title.y = element_blank(),
              strip.background = element_blank(),
              strip.text.y = element_blank()) + 
        # draw_image(here("material", "20180712_glyc_glycTMG20uM_Pos0_t151.jpg"), scale=2, x=0.4, y=1) +
        # draw_image(here("material", "20181008_glyc_lactuloseTMG20uM_Pos6_t151.jpg"), scale=2) +
        NULL, 
      nrow=1, labels=c('F', 'G', '', ''), rel_widths=c(0.7, 0.3, 0.85, 1.15)),
    ncol=1, rel_heights=c(1.3, 0.1, 1.0))
}) ()

save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]](),
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.2
)


####

# (myfigs[[4]] <- function(){
#   pdftools_installed <- require(pdftools)
#   plot_grid(
#     if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1), scale=1.0),
#     myplots[['autoactiv_model_phase']] +
#       # theme(axis.title.y = element_text(margin=margin(l=-2, r=4))) +
#       NULL,
#     if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3), scale=1.0),
#     myplots[['twocmps_model_phase']],
#     labels="AUTO", ncol=1, rel_heights=c(1, 1, 1), align='v'
#   )}
# )()
# 
# save_plot(here("plots", "MoM_lacDilution_fig4.pdf"), myfigs[[4]](),
#           base_height=NULL, base_width=2.25 * 14/8, # 2 cols
#           base_aspect_ratio = 1/2
# )

# (myfigs[[4]] <- function(){
#   pdftools_installed <- require(pdftools)
#   plot_grid(
#     ggdraw() +
#       ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1), 
#                                                        x=-.105, y=-.1, scale=.4) ) +
#       draw_plot(myplots[['autoactiv_model_phase']]) +
#       NULL,
#     ggdraw() +
#       ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3), 
#                                                       x=-.105, y=-.1, scale=.4) ) +
#       draw_plot(myplots[['twocmps_model_phase']]) +
#       NULL,
#     labels="AUTO", ncol=1, rel_heights=c(1, 1), align='v'
#   )}
# )()

(myfigs[[4]] <- function(){
  pdftools_installed <- require(pdftools)
  plot_grid(
    NULL, 
    myplots[['autoactiv_model_induction']] +
      scale_y_log10(breaks=c(1e-3, 1e0, 1e3),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      expand_limits(y=1.e-4) +
      # expand_limits(y=2e-3) +
      # theme(axis.title.y = element_text(margin=margin(l=-2, r=4))) +
      ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1),
                                                      x=1.5, y=-2.8, scale=4.35) ) +
      NULL,
    myplots[['autoactiv_model_phase']],
    NULL,
    myplots[['twocmps_model_phase']] +
      ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3),
                                                      x=-.2, y=-3.7, scale=5) ) +
      NULL,
    labels=c("A  auto-activating operon", "", "B", "C  two-components system", ""), 
    ncol=1, rel_heights=c(.1, 1, 1, .1, 1.2), align='v', hjust=-.05
  )}
)()

save_plot(here("plots", "MoM_lacDilution_fig4.pdf"), myfigs[[4]](),
          base_height=NULL, base_width=2.25 * 14/8, # 2 cols
          base_aspect_ratio = 1/2
)

