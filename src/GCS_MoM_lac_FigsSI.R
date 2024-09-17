
mytables[['lacl_list']] %>% 
  mutate(n_lags=ifelse(condition %in% c("switch_glycerol_TMG20", "switch_lactulose"), NA, n_lags)) %>% 
  ungroup() %>% select(-condition, -n_div_cells, -n_obs) %>% rename(condition=label) %>%
  mutate(condition=str_replace(condition, '>', ' to ')) %>%
  (function(.df)
    knitr::kable(.df, "latex", booktabs=TRUE, #longtable = TRUE,
                 label="lactulose-list", 
                 col.names=c('condition', 'date', '# growth channels', '# cells at switch', '# estimated lags', '# arrested cells at switch'),
                 caption='List of experiments on growth-coupled sensitivity during transient growth arrest (Fig. 2D) with summary statistics.') %>%
     # kableExtra::kable_styling(full_width=TRUE) %>%
     kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
     kableExtra::column_spec(3:6, width = "3.5em") %>%
     kableExtra::column_spec(7:8, width = "5em") %>%
     kableExtra::row_spec(which(.df$date %in% discarded_dates), italic=T, color="gray") %>%
     identity()
  ) %>%
  # str_replace(fixed("{tab:}"), "{tab:lactulose-list}") %>%
#   str_replace(fixed("\\resizebox{\\linewidth}{!}"), "\\resizebox*{!}{0.9\\textheight}") %>% 
  write(here('plots', 'SI_figs', 'lacl-list.tex'))


(myplots[['TMG_induction_gly04_hill']] <- (function() {
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
    scale_y_log10() +
    labs(x="TMG concentration (µM)", y="[LacZ-GFP] per cell (AU)    ") +
    NULL
  })())

(myplots[['TMG_induction_gly04']] <- (function() {
  load('data/20180703_ASC662_M9gly04pc_TMG.RData')
  mydata %>% ungroup() %>% 
    filter(Ch==2) %>% 
    ggplot(aes(tmg, gfp)) +
    geom_vline(xintercept = 20, lty="dashed", col=qual_cols[2], size=0.8) + 
    geom_point(alpha=.2, stroke=0, position=position_jitter(width=.1)) +
    scale_x_continuous(trans='log2', limits=c(5, 80), breaks=c(6, 12, 25, 50)) +
    scale_y_log10() +
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
  ncol=1, labels=c("", "C"), rel_heights=c(1, 1.2)
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
# MILLER ASSAY ####
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

plot_grid(
  myplots[['Cline_You2013']] +
    scale_x_continuous(breaks=c(0, 0.5, 1)) +
    coord_cartesian(xlim=c(0, 1.25), ylim=c(0, NA)) +
    expand_limits(x=c(-.1, 1.4)) +
    labs(y='max [LacZ] (MU)      ') +
    NULL,
  myplots[['lac_model_crp_phdiag']] +
    labs(y="external [TMG] (µM)"),
  # NULL,
  nrow=1, labels="AUTO") %>% 
  save_plot(here("plots", "SI_figs", "native-lac-model.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 2.7)
  


###### #
# SUGAR MIXTURES ####

mytables[['sugarmix_list']] %>%
  mutate(treatment=fct_recode(treatment, 'glucose only'='none', 'with lactose 0.55 mM'='lac002', 'with IPTG 200µM'='iptg')) %>% 
  mutate_if(is.numeric, list(as.character)) %>% # protect existing rounding of floats
  knitr::kable("latex", booktabs=TRUE, #longtable = TRUE,
               label="sugarmix-list", 
               col.names=c('treatment', '[glucose] (µM)', 'date', '# analysed cells', 'prop. induced', 
                           '# growth channels', '# frames'),
               caption='List of experiments on the concentration-dependent sugar preference (Fig. 3B-C) with summary statistics. Experiments that have been discarded from further analysis are greyed out.') %>%
  # kableExtra::kable_styling(full_width=TRUE) %>%
  kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
  # kableExtra::column_spec(3:6, width = "3.5em") %>%
  # kableExtra::column_spec(7:8, width = "5em") %>%
  kableExtra::row_spec(
    which(mytables[['sugarmix_list']]$date==20210122 & mytables[['sugarmix_list']]$glu==16
          | mytables[['sugarmix_list']]$date %in% c(20210122, 20210305) & mytables[['sugarmix_list']]$glu==32
    ), italic=T, color="gray") %>%
  # str_replace(fixed("{tab:}"), "{tab:sugarmix-list}") %>%
  write(here('plots', 'SI_figs', 'sugarmix-list.tex'))


plot_grid(
  plot_grid(
    myplots[['sugarmix_m9zero']] +
      labs(title="no C source") ,
    myplots[['sugarmix_indn_thr']] +
      labs(col="[glucose]\n(µM)") +
      theme_cowplot_legend_inset() +
      labs(title="with lactose 0.55 mM") ,
  nrow=1, rel_widths = c(1, 1.2), labels=c("A", "C")),
  myplots[['sugarmix_gr_pos']] +
    # guides(col=guide_colorbar()) +
    labs(title="glucose only", col="cell rank in channel") +
    theme(legend.position = "top"),
  ncol=1, rel_heights = c(1, 3), labels=c("", "B")) %>% 
  save_plot(here("plots", "SI_figs", "sugarmix-controls.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 1 col
            base_aspect_ratio = 1)

plot_grid(
  myplots[['sugarmix_monod']] +
    scale_x_continuous(breaks=seq(0, 2000, by=500)),
  myplots[['sugarmix_monod_inv']],
  nrow=1, rel_widths = c(1, 1), labels="AUTO") %>% 
  save_plot(here("plots", "SI_figs", "sugarmix-monod.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 1 col
            base_aspect_ratio = 2.5)

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



#### XXXX TO CLEAN UP

cowplot::plot_grid(
  myplots[['sugarmix_singlecell_rates']], 
  myplots[['sugarmix_singlecell_cdfs']],
  labels = c("A","B"), nrow = 2, rel_heights = c(0.55, 0.45)
) %>% 
save_plot(here("plots", "SI_figs", "sugarmix-switching.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/7, # 2 cols
          base_aspect_ratio = 1.4)


cowplot::plot_grid(
  nrow = 2, rel_heights = c(0.55, 0.45), labels = c("A",""),
  myplots[['sugarmix_singlecell_rates']], 
  cowplot::plot_grid(
    ncol = 2, rel_widths = c(0.6, 0.4), labels = c("B","C"), 
    myplots[['sugarmix_singlecell_cdfs']],
    myplots[['sugarmix_singlecell_pval']])
  ) %>% 
  save_plot(here("plots", "SI_figs", "sugarmix-switching-3panels.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/7, # 2 cols
            base_aspect_ratio = 1.4)

