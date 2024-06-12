library(mvtnorm)
library(ggplot2)

source("/accelerated-approval/functions.R")

# Figure 1 of the futility boundary for the success boundary 0.01, with different correlation values
# and information fractions
alpha <- 0.025
p_S <- 0.01
r <- seq(0.1, 0.7, 0.2)
rho <- seq(-1, 1, 0.001)
scen <- expand.grid(rho = rho, r = r)
p_F <- rep(NA, nrow(scen))
for (i in 1:nrow(scen)) {
  p_F[i] <- uniroot(type1_fun_alpha,
                   c(p_S, 1),
                   alpha = alpha,
                   r = scen$r[i],
                   rho = scen$rho[i],
                   p_S = p_S)$root
}
data <- data.frame(scen, p_F)
data$r <- factor(data$r)

ggplot(data, aes(x = rho, y = p_F, group = r)) +
  geom_line(aes(linetype = r), linewidth = 1.2) +
  scale_x_continuous(breaks = seq(-1, 1, 0.2),
                     limits = c(-1, 1),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(seq(0, 1, 0.1), alpha),
                     minor_breaks = c(seq(0, 1, 0.1), alpha),
                     limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_linetype_manual(values=c("solid", "longdash", "dotdash", "dotted"))+
  xlab("Correlation between ORR and OS") +
  ylab("Futility boundary p_F for ORR") +
  guides(linetype = guide_legend(title = "Information\nfraction of OS",
                                 position = "inside")) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position.inside = c(0.85, 0.79),
        legend.key.width = unit(2, 'cm'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13))

# Table 1 of the futility boundary for the success boundary 0.01
alpha <- 0.025
p_S <- 0.01
r <- seq(0.1, 0.7, 0.1)
p_F <- r
orr_c <- 0.2
orr_t <- r
nper <- 100
orr_t_0.025 <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = alpha)$root
for (i in 1:length(r)) {
  p_F[i] <- uniroot(type1_search, c(p_S, 1), alpha = alpha, r = r[i], p_S = p_S)$root
  orr_t[i] <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = p_F[i])$root
}
rbind(r = formatC(r, format = "f", digits = 1),
      p_F = formatC(p_F, format = "f", digits = 4),
      diff = formatC(round((orr_t - orr_c) / (orr_t_0.025 - orr_c) *100, 0), format = "f", digits = 0)
)

# Table 2 of the futility boundary for the success boundary 0.005
alpha <- 0.025
p_S <- 0.005
r <- seq(0.1, 0.7, 0.1)
p_F <- r
orr_c <- 0.2
orr_t <- r
nper <- 100
orr_t_0.025 <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = alpha)$root
for (i in 1:length(r)) {
  p_F[i] <- uniroot(type1_search, c(p_S, 1), alpha = alpha, r = r[i], p_S = p_S)$root
  orr_t[i] <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = p_F[i])$root
}
rbind(r = formatC(r, format = "f", digits = 1),
      p_F = formatC(p_F, format = "f", digits = 4),
      diff = formatC(round((orr_t - orr_c) / (orr_t_0.025 - orr_c) *100, 0), format = "f", digits = 0)
)

# Table 2 of the futility boundary for the success boundary 0.0125
alpha <- 0.025
p_S <- 0.0125
r <- seq(0.1, 0.7, 0.1)
p_F <- r
orr_c <- 0.2
orr_t <- r
nper <- 100
orr_t_0.025 <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = alpha)$root
for (i in 1:length(r)) {
  p_F[i] <- uniroot(type1_search, c(p_S, 1), alpha = alpha, r = r[i], p_S = p_S)$root
  orr_t[i] <- uniroot(orr_fun, c(0, 1), orr_c = orr_c, nper = nper, alpha = p_F[i])$root
}
rbind(r = formatC(r, format = "f", digits = 1),
      p_F = formatC(p_F, format = "f", digits = 4),
      diff = formatC(round((orr_t - orr_c) / (orr_t_0.025 - orr_c) *100, 0), format = "f", digits = 0)
)


