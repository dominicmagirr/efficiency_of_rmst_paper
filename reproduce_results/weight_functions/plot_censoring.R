t_R <- 0.5
t_H <- 3
t_seq <- seq(0, t_H, by = 0.01)
y <- ifelse(t_seq < t_H - t_R, 1, (t_H - t_seq) / t_R)

plot(t_seq, y, type = "l")

## ggplot

## prob (censoring)