
mytables[['lacl_list']] %>% 
  mutate(n_lags=ifelse(condition %in% c("switch_glycerol_TMG20", "switch_lactulose"), NA, n_lags)) %>% 
  ungroup() %>% select(-condition) %>% rename(condition=label) %>%
  mutate(condition=str_replace(condition, '>', ' to ')) %>%
  (function(.df)
    knitr::kable(.df, "latex", booktabs=TRUE, #longtable = TRUE,
                 col.names=c('condition', 'date', '# growth channels', '# full cell cycles', '# observations',
                             '# cells at switch', '# estimated lags', '# arrested cells at switch'),
                 caption='List of experiments on growth-coupled sensitivity during transient growth arrest (Fig. 1C) with summary statistics.') %>%
     # kableExtra::kable_styling(full_width=TRUE) %>%
     kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
     kableExtra::column_spec(3:6, width = "3.5em") %>%
     kableExtra::column_spec(7:8, width = "5em") %>%
     kableExtra::row_spec(which(.df$date %in% discarded_dates), italic=T, color="gray") %>%
     identity()
  ) %>%
  str_replace(fixed("{tab:}"), "{tab:lactulose-list}") %>%
#   str_replace(fixed("\\resizebox{\\linewidth}{!}"), "\\resizebox*{!}{0.9\\textheight}") %>% 
  write(here('plots', 'SI_figs', 'lacl-list.tex'))


(myplots[['TMG_induction_gly04']] <- (function() {
  load('data/20180703_ASC662_M9gly04pc_TMG.RData')
  mydata %>% 
    ungroup() %>% 
    filter(Ch==2) %>% 
    ggplot(aes(tmg, gfp)) +
    geom_point(alpha=.2, stroke=0, position=position_jitter(width=.1)) +
    stat_function(fun=hill_fn, args=as.list(summary(myfit_all)$coefficients[,1]), col="red") +
    stat_function(fun=hill_fn, args=as.list(summary(myfit_all)$coefficients[,1]), col="red", 
                  geom="point", n=1, xlim=log2(c(20, 20)), size=3) +
    scale_x_continuous(trans='log2', limits=c(5, 200), breaks=c(6, 12, 25, 50, 100, 200)) +
    scale_y_continuous(trans='log2', breaks=c(16, 128, 1024)) +
    labs(x="TMG concentration (µM)", y="[LacZ-GFP] per cell (AU)    ") +
    NULL
  })())
  # save_plot(here("plots", "SI_figs", "induction-tmg-gly.pdf"), .,
  #           base_height=NULL, base_width=4.75 * 14/7, # 2 cols
  #           base_aspect_ratio = 2)

# plot_grid(
#   myplots[['TMG_induction_gly04']],
#   myplots[['lags_hist_lacl']],
#   nrow=1, labels=c('A', 'B'), align="h"
# )
plot_grid(
  plot_grid(
    myplots[['TMG_induction_gly04']],
    myplots[['lags_hist_lacl']],
    nrow=1, labels=c('A', 'B'), align="h"
  ),
  myplots[['TMG_switch_gr_hist']] +
    theme_cowplot_legend_inset(),
  ncol=1, labels=c("", "B"), rel_heights=c(1, 1.2)
) %>% 
  save_plot(here("plots", "SI_figs", "transient-arrest-tmg.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.6)

plot_grid(
  plot_grid(
    myplots[['lacl_gfp_facets']] + theme(legend.position = 'top'),
    myplots[['lacl_gr_facets']] + theme(legend.position = 'top'),
    ncol=1, rel_heights = c(5, 4), labels=c('A', 'B'), align='v'),
  myplots[['lacl_gr_lacz']] +
    facet_wrap(~condition, ncol=1, labeller=as_labeller(rename_conds)),
  nrow=1, rel_widths = c(1.6, 1), labels=c("", "C")
) %>% 
  save_plot(here("plots", "SI_figs", "lacl-gr-lacz.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1)



############
# MILLER ASSAY
plot_grid(
  myplots_miller[['GCS_MA_examples']] + 
    theme_cowplot_legend_inset(),
  myplots_miller[['GCS_calib']] +
    theme_cowplot_legend_inset(),
  ncol=1, rel_heights = c(1, 1.7), labels="AUTO"
) %>% 
  save_plot(here("plots", "SI_figs", "miller-example-calib.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.2)

(myplots_miller[['GCS_sugars_gr']] +
    guides(col=guide_legend(nrow=2, byrow=TRUE)) +
    theme_cowplot_legend_inset() +
  NULL) %>% 
  save_plot(here("plots", "SI_figs", "miller-sugars-gr.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.4)

(myplots_miller[['GCS_sugars_indct']]() +
    guides(col=guide_legend(nrow=2, byrow=TRUE)) +
    theme_cowplot_legend_inset() +
    NULL) %>% 
  save_plot(here("plots", "SI_figs", "miller-sugars-indn.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.4)

(myplots_miller[['GCS_sugars']]() +
    theme_cowplot_legend_inset() +
    guides(col=guide_legend(ncol=2, title.position = "top")) +
    NULL) %>% 
  save_plot(here("plots", "SI_figs", "miller-sugars-gcs.pdf"), .,
            base_height=NULL, base_width=2.25 * 14/7, # 1 cols
            base_aspect_ratio = 0.9)

(myplots[['GCS_cm_gr']] <- function() {
  load(list.files(c("share", "data"), "GCS_Miller_lac_cm_\\w+\\.RData", full.names=TRUE) ) # do NOT use load() after %>%
  mygcs %>%
    filter(!is.na(replicate)) %>%
    filter(elapsed>-2*3600, od>5e-4) %>% 
    left_join(myplates %>% select(date, strain, media, row, tmg) %>% distinct()) %>% 
    filter(strain=='MG1655') %>%
    group_by(date, strain, media, tmg, replicate) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(mod_gc=map(data, ~lm(log(od)~elapsed, data=.)),
           gr=map_dbl(mod_gc, ~coef(.)[2]), gr_se=map_dbl(mod_gc, ~summary(.)$coefficients[2,2]),
           od_log=map_dbl(mod_gc, ~predict(., data.frame(elapsed=5*60)) ), od=exp(od_log),
           date=paste('day', as.numeric(factor(date))),
    ) %>% 
    filter(gr*3600/log(2) <  1.4) %>% 
    ggplot(aes(tmg, gr*3600/log(2), col=factor(date))) +
    facet_wrap(~media) + 
    geom_hline(yintercept = 0, lty='dotted') +
    # geom_vline(xintercept = 50) +
    geom_pointrange(aes(ymin=gr-gr_se, ymax=gr+gr_se, shape=replicate), position=position_jitter(width=.025)) +
    scale_x_log10(limits=c(3, NA), breaks=c(10, 100, 1000)) + 
    ylim(-0.1, NA) +
    # scale_color_discrete() +
    labs(x="[TMG] (µM)", y="growth rate (dbl/h)", col="date") +
    guides(col=guide_legend(nrow=1, byrow=TRUE)) +
    theme_cowplot_legend_inset() +
    theme(legend.position = 'top') +
    NULL
} )() %>% 
  save_plot(here("plots", "SI_figs", "miller-cm-gr.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.7)

(myplots[['GCS_cm_indct']] <- function() {
  load(list.files(c("share", "data"), "GCS_Miller_lac_cm_\\w+\\.RData", full.names=TRUE) ) # do NOT use load() after %>%
  myindct %>% unnest(data) %>%
    filter(strain=='MG1655') %>% 
    ungroup() %>% 
    mutate(media=factor(media, levels=sort(levels(media))),
           date=paste('day', as.numeric(factor(date))),
    ) %>% 
    # myindct %>% unnest(data_fit) %>%
    ggplot(aes(tmg, activity, col=date)) +
    facet_wrap(~media) +
    geom_point(aes(shape=replicate)) +
    geom_line(data=unnest(myindct, prediction) %>% filter(strain=='MG1655') %>%  ungroup() %>% 
                mutate(date=paste('day', as.numeric(factor(date)))) ) +
    scale_x_log10(limits=c(3, NA), breaks=c(1, 10, 100, 1000)) +
    scale_y_log10(breaks=c(1, 30, 1000)) +
    labs(x='TMG concentration (µM)', y='p_lac activity (AU)') +
    guides(col=guide_legend(nrow=1, byrow=TRUE)) +
    theme_cowplot_legend_inset() +
    theme(legend.position = 'top') +
    NULL
} )() %>% 
  save_plot(here("plots", "SI_figs", "miller-cm-indn.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.7)

# (myplots_miller[['GCS_cm_maxindct']](strain=='MG1655') +
#     guides(shape='none') +
#     NULL) %>% 
#   save_plot(here("plots", "SI_figs", "miller-cm-maxindn.pdf"), .,
#             base_height=NULL, base_width=2.25 * 14/7, # 1 col
#             base_aspect_ratio = 1)

###### #
# SUGAR MIXTURES ####

mytables[['sugarmix_list']] %>%
  mutate(treatment=fct_recode(treatment, 'glucose only'='none', 'with lactose 200mg/L'='lac002', 'with IPTG 200µM'='iptg')) %>% 
  mutate_if(is.numeric, list(as.character)) %>% # protect existing rounding of floats
  knitr::kable("latex", booktabs=TRUE, #longtable = TRUE,
               col.names=c('treatment', '[glucose] (mg/L)', 'date', '# analysed cells', 'prop. induced', 
                           '# growth channels', '# frames'),
               caption='List of experiments on the concentration-dependent sugar preference (Fig. 3) with summary statistics. Experiments that have been discarded from further analysis are greyed out.') %>%
  # kableExtra::kable_styling(full_width=TRUE) %>%
  kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
  # kableExtra::column_spec(3:6, width = "3.5em") %>%
  # kableExtra::column_spec(7:8, width = "5em") %>%
  kableExtra::row_spec(
    which(mytables[['sugarmix_list']]$date==20210122 & mytables[['sugarmix_list']]$glu==2.9
          | mytables[['sugarmix_list']]$date %in% c(20210122, 20210305) & mytables[['sugarmix_list']]$glu==5.8
    ), italic=T, color="gray") %>%
  str_replace(fixed("{tab:}"), "{tab:sugarmix-list}") %>%
  write(here('plots', 'SI_figs', 'sugarmix-list.tex'))


plot_grid(
  plot_grid(
    myplots[['sugarmix_m9zero']] +
      labs(title="no C source") ,
    myplots[['sugarmix_indn_thr']] +
      labs(col="[glucose]\n(mg/L)") +
      theme_cowplot_legend_inset() +
      labs(title="with lactose 200 mg/L") ,
  nrow=1, rel_widths = c(1, 1.2), labels=c("A", "C")),
  myplots[['sugarmix_gr_pos']] +
    # guides(col=guide_colorbar()) +
    labs(title="glucose only", col="cell rank in channel") +
    theme(legend.position = "top"),
  ncol=1, rel_heights = c(1, 3), labels=c("", "B")) %>% 
  save_plot(here("plots", "SI_figs", "sugarmix-controls.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 1 col
            base_aspect_ratio = 1)

(myplots[['sugarmix_gr_distr']] +
    theme_cowplot_legend_inset() +
    NULL) %>% 
  save_plot(here("plots", "SI_figs", "sugarmix-gr-distr.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 1 col
            base_aspect_ratio = 0.9)

# myplots[['sugarmix_crp']] %>% 
#   save_plot(here("plots", "SI_figs", "sugarmix-crp.pdf"), .,
#             base_height=NULL, base_width=2.25 * 14/7, # 1 col
#             base_aspect_ratio = 1)
