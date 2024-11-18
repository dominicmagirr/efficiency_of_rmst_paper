if (!dir.exists("figs")) {
  dir.create("figs")
}

cowplot::ggsave2("figs/censoring.pdf", 
       cowplot::plot_grid(p_1, p_2, p_3, p_4), 
       width = 9, 
       height = 6, 
       units = "in", 
       dpi = 300)
