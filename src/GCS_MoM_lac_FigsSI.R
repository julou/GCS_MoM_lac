
#  xxx add comments + flow control
# remove lac_iptg
mytables[['expts_list']] %>% 
  filter(! condition %in% c("switch_ramp15min", "switch_lactose_priming"),
         ! str_detect(condition, "_stdIllum")) %>% 
  ungroup() %>% select(-condition) %>% rename(condition=label) %>% 
  mutate(condition=str_replace(condition, '>', ' to ')) %>% 
  (function(.df)
    knitr::kable(.df, "latex", booktabs=TRUE, #longtable = TRUE, 
                 col.names=c('condition', 'date', '# growth channels', '# full cell cycles', '# observations',
                             '# cells at switch', '# estimated lags', '# arrested cells at switch'),
                 caption='List of experiments used in this study with summary statistics. Experiments that have been discarded from further analysis are greyed out.') %>% 
     # kableExtra::kable_styling(full_width=TRUE) %>% 
     kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
     kableExtra::column_spec(3:6, width = "3.5em") %>% 
     kableExtra::column_spec(7:8, width = "5em") %>% 
     kableExtra::row_spec(which(.df$date %in% discarded_dates), italic=T, color="gray") %>% 
     identity()
  ) %>% 
  str_replace(fixed("{tab:}"), "{tab:expts-list}") %>% 
  str_replace(fixed("\\resizebox{\\linewidth}{!}"), "\\resizebox*{!}{0.9\\textheight}") %>% 
  write(here('plots', 'SI_figs', 'expts-list.tex'))



myplots[['gr_length_medians']](.article_ds=TRUE) %>% 
  save_plot(here("plots", "SI_figs", "gr-length-medians.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.35)

# plot_grid(
#   myplots[['lags_gr_before']](.size_n=FALSE) +
#     theme(legend.position=c(1,1), legend.justification=c(1,1),
#           axis.title.x=element_blank()) +
#     NULL,
#   myplots[['gr_hist_before']] +
#     theme(legend.position='none') +
#     NULL,
#   ncol=1, rel_heights=c(1.5, 1), align='v'
# )
# # gr-before-filtering


myplots[['naive_lags_per_expt']]() %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-expt.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.2)


myplots[['naive_lags_per_pos']] %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-pos.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.6)

myplots[['lags_hist_ramp']] %>% 
  save_plot(here("plots", "SI_figs", "lags-hist-ramp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)

(myplots[['glyc_mix_violin']] +
  theme(legend.position = 'none')) %>% 
  save_plot(here("plots", "SI_figs", "glyc-mix-violin.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2)

myplots[['basal_perturb_gfp']] %>% 
  save_plot(here("plots", "SI_figs", "basal-perturb-gfp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)

myplots[['naive_lags_correl']] %>% 
save_plot(here("plots", "SI_figs", "naive-lags-correl.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2)

# (myplots[['high-illum-sensitivity']] <-
# bind_rows(
#   read.csv("../data/20180403_glu_lac_hiIllum_1_pos2_t215_MG1655_hist.csv") %>%
#     with(rep.int(value, count)) %>% tibble(strain='MG1655', int=.),
#   read.csv("../data/20180403_glu_lac_hiIllum_1_pos2_t215_ASC662_hist.csv") %>%
#     with(rep.int(value, count)) %>% tibble(strain='ASC662', int=.)
#   ) %>%
#   ggplot(aes(int, ..density.., col=strain)) +
#   geom_freqpoly(binwidth=50) +
#   scale_colour_discrete(breaks=c('MG1655', 'ASC662'), labels=c('MG1655 (70 cells)', 'ASC662 (93 cells)')) +
#   labs(x='pixel intensity (AU)') +
#   NULL)

plot_grid(
  plot_grid(
    myplots[['lags_gfp_scatter']](.fit=FALSE) +
      theme(legend.position = 'none',
            axis.title.x = element_blank()),
    myplots[['lags_gfp_diff_cdf']] +
      scale_y_continuous(breaks=c(0, 0.5, 1)) +
      theme(legend.position = 'none'),
    ncol=1, rel_heights = c(1.2, 0.8), align='v'),
  get_legend(myplots[['lags_gfp_diff_cdf']]) ,
  nrow=1, rel_widths = c(0.85, 0.15)
) %>% 
  save_plot(here("plots", "SI_figs", "lag-gfp-ini.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.35)



# myplots[['memory_cdfs_facets']] %>% 
#   save_plot(here("plots", "SI_figs", "memory-cdfs-facets.pdf"), .,
#             base_height=NULL, base_width=4.75 * 14/8, # 2 cols
#             base_aspect_ratio = 1.35)


myplots[['lag_memory']] %>% 
  save_plot(here("plots", "SI_figs", "lag-memory.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.2)


myplots[['memory_elapsed_divs']] %>% 
  save_plot(here("plots", "SI_figs", "memory-elapsed-divs.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


myplots[['lags_inherited_gfp']] %>% 
  save_plot(here("plots", "SI_figs", "lags-inherited-gfp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.3)

plot_grid(
  myplots[['naive_arrest_cdf']] +
     annotate('segment', x=15, xend=15, y=0.99, yend=0.01,
              arrow=arrow(length=unit(0.1, "inches")), lineend='butt', linejoin='mitre') +
     annotate('text', x=15, y=-.05, label="15", size=4) +
     coord_cartesian(ylim=c(0, 1.02), clip='off') +
     NULL,
   myplots[['lags_types_correl']],
   nrow=1) %>% 
save_plot(here("plots", "SI_figs", "growth-lags.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2.5)

(myplots[['TMG_induction_gly04']] <- (function() {
  load('data/20180703_ASC662_M9gly04pc_TMG.RData')
  mydata %>% 
    filter(Ch==2) %>% 
    ggplot(aes(tmg, gfp)) +
    geom_point(alpha=.2, stroke=0, position=position_jitter(width=.1)) +
    stat_function(fun=function(c, ...) log(hill_fn(c, ...)) / log(2), args=as.list(summary(myfit_all)$coefficients[,1]), col="red") +
    stat_function(fun=function(c, ...) log(hill_fn(c, ...)) / log(2), args=as.list(summary(myfit_all)$coefficients[,1]), col="red", 
                  geom="point", n=1, xlim=log2(c(20, 20)), size=3) +
    scale_x_continuous(trans='log2', limits=c(5, 200), breaks=c(6, 12, 25, 50, 100, 200)) +
    scale_y_continuous(trans='log2', breaks=c(16, 128, 1024)) +
    labs(x="TMG concentration (µM)", y="LacZ-GFP concentration per cell (AU)") +
    NULL
  })()) %>% 
  save_plot(here("plots", "SI_figs", "induction-tmg-gly.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


myplots[['TMG_switch_gr_hist']] %>% 
  save_plot(here("plots", "SI_figs", "growth-arrest-tmg.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


plot_grid(
  myplots[['lacl_gfp_facets']] + theme(legend.position = 'right'),
  myplots[['lacl_gr_facets']] + theme(legend.position = 'right'),
  ncol=1, rel_heights = c(5, 4), align='v') %>% 
  save_plot(here("plots", "SI_figs", "lacl-gfp-gr-hists.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1)


myplots[['lacl_gr_lacz']] %>% 
  save_plot(here("plots", "SI_figs", "lacl-gr-lacz.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)

(myplots[['signal_dep_decay']] <-
    tibble(path=list.files(here("material"), "signal_coupled_decay_pd_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest() %>% 
    extract(path, c('curve'), ".*/signal_coupled_decay_pd_(.+)\\.txt") %>% 
    # with(range(X1))
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                   tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) ))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    # geom_hline(yintercept = 0.014269, lty='dashed') +
    # annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(limits=c(.5, 1e3), breaks=c(1, 24, 500), expand=c(0, 0)) +
    scale_y_log10(limits=c(8, 1e5), breaks=c(1e2, 1e4), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.35) +
    # labs(x='doubling time (h)', y=expression(paste('signal strength ', italic('s'), '/', italic('s[0]')))) +
    labs(x='doubling time (h)', y=expression(paste('signal strength ', s / s[0] ))) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

# (function() {
#   pdftools_installed <- require(pdftools)
#   myplots[['signal_dep_decay']] +
#     # theme(axis.title.y = element_text(margin=margin(l=-2, r=4))) +
#     ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=2),
#                                                     x=-.1, y=.8, scale=2.5) ) +
#     NULL
#   
# })() %>%
#   save_plot(here("plots", "SI_figs", "signal-dep-decay.pdf"), .,
#             base_height=NULL, base_width=2.25 * 14/8, # 1 col
#             base_aspect_ratio = 1/0.75
#   )

(function() {
  pdftools_installed <- require(pdftools)
  plot_grid(
    if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=2) ),
    myplots[['signal_dep_decay']],
    nrow=1)  
})() %>% 
  save_plot(here("plots", "signal-dep-decay.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.5
  )


