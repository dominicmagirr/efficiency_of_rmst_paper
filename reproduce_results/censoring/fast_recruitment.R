#########################################################################
## recruitment assumptions

dat <- data.frame(n_pats = c(0, 30, 60, 90, 120, 150, 150),
                  cal_t = c(0, .1, .2, .3, .4, .5, 3))

library(ggplot2)

dat_lines <- data.frame(x = c((0:4)/10, rep(3, 5)),
                        id = rep(as.character(0:4), 2),
                        y = rep((0:4)*30, 2)) 

dat_ribbon <- data.frame(x = c(0, 0.45),
                         ymin = c(10, 150),
                         ymax = c(150, 150),
                         n_pats = c(0, 0.45))

p_1 <- ggplot(data = dat,
       mapping = aes(x = cal_t,
                     y = n_pats)) +
  geom_line() + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(data = dat_lines,
            mapping = aes(x = x, 
                          y = y,
                          group = id),
            arrow = arrow(),
            linetype = 2) +
  xlab("Calendar time") +
  ylab("Number of patients recruited") + 
  ggtitle("Fast recruitment") + 
  geom_ribbon(data = dat_ribbon,
              mapping = aes(x = x,
                            ymin = ymin,
                            ymax = ymax),
              fill = "grey") +
  geom_vline(xintercept = 3)


###########################################################################
## censoring distribution


dat <- data.frame(pat_t = c(0, 2.5, 3),
                  p_cens = c(1, 1, 0))

dat_ribbon <- data.frame(x = c(2.55, 3),
                         ymin = c(1, 0.1),
                         ymax = c(1, 1))



p_2 <- ggplot() +
  geom_line(data = dat,
            mapping = aes(x = pat_t,
                          y = p_cens)) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Patient time") +
  ylab("1 - P(Censored)") + 
  ggtitle("Light censoring") + 
  geom_ribbon(data = dat_ribbon,
              mapping = aes(x = x,
                            ymin = ymin,
                            ymax = ymax),
              fill = "grey") 




























