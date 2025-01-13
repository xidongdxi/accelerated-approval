library(mvtnorm)
library(ggplot2)
library(latex2exp)

source("/accelerated-approval/functions.R")

# Figure 1 of the futility boundary for the success boundary 0.01, with different correlation values
# and information fractions
alpha <- 0.025
p_S <- 0.005
r <- seq(0.1, 0.9, 0.2)
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

pdf(file = "C:\\Users\\dxi1\\OneDrive - Gilead Sciences\\Paper\\AA\\pF_rho.pdf",
    width = 8,
    height = 5)

ggplot(data, aes(x = rho, y = p_F, group = r)) +
  geom_line(aes(linetype = r), linewidth = 1.2) +
  scale_x_continuous(breaks = seq(-1, 1, 0.2),
                     limits = c(-1, 1),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = c(seq(0, 1, 0.1), alpha),
                     minor_breaks = c(seq(0, 1, 0.1), alpha),
                     limits = c(0, 1),
                     expand = c(0.01, 0.01)) +
  scale_linetype_manual(values=c("solid", "longdash", "dashed", "dotdash", "dotted"))+
  xlab(TeX(r'(Correlation  $\rho$  between test statistics of ORR and OS)')) +
  ylab(TeX(r'(Futility $p$-value boundary  $p_F$  for ORR)')) +
  guides(linetype = guide_legend(title = paste("Information\nfraction ", TeX(r'($r$)'), " of OS"),
                                 position = "inside")) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position.inside = c(0.15, 0.3),
        legend.key.width = unit(2, 'cm'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13))

dev.off()

# Table 1 of the futility boundary for the success boundary 0.005
alpha <- 0.025
p_S <- 0.005
r <- seq(0.1, 0.9, 0.1)
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
      p_F = formatC(p_F, format = "f", digits = 3),
      diff = paste0(formatC(round((orr_t - orr_c) *100, 1), format = "f", digits = 1), "%"),
      diff_rel = paste0(formatC(round((orr_t - orr_c) / (orr_t_0.025 - orr_c) *100, 0), format = "f", digits = 0), "%")
)

# Table 2 of the futility boundary for the success boundary 0.01
alpha <- 0.025
p_S <- 0.01
r <- seq(0.1, 0.9, 0.1)
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
      p_F = formatC(p_F, format = "f", digits = 3),
      diff = paste0(formatC(round((orr_t - orr_c) *100, 1), format = "f", digits = 1), "%"),
      diff_rel = paste0(formatC(round((orr_t - orr_c) / (orr_t_0.025 - orr_c) *100, 0), format = "f", digits = 0), "%")
)

###############################################
# Impact on OS power from the futility boundary
library(mvtnorm)
alpha <- 0.025
p_S <- 0.005
alpha_ORR <- p_S
alpha_OS <- alpha - alpha_ORR
p_F <- 0.254
r <- 0.5
power_ORR <- 0.965
power_OS <- 0.9

ncp_ORR <- qnorm(1 - alpha) - qnorm(1 - power_ORR)
ncp_OS <- qnorm(1 - alpha) - qnorm(1 - power_OS)

pnorm(qnorm(1 - alpha_ORR), ncp_ORR, 1, lower.tail = F)

loss <- NULL
power_OS <- NULL
for (rho in seq(-1, 1, 0.001)) {
  cr <- matrix(c(1, rho * sqrt(r), rho * sqrt(r), 1), nrow = 2)
  # ORR not failed and OS successful
  power_OS <- c(power_OS, pmvnorm(upper = c(Inf, Inf),
                          lower = c(qnorm(1 - p_F), qnorm(1 - alpha)),
                          mean = c(ncp_ORR, ncp_OS),
                          corr = cr))
  # ORR failed but OS successful
  loss <- c(loss, pmvnorm(upper = c(qnorm(1 - p_F), Inf),
                          lower = c(-Inf, qnorm(1 - alpha)),
                          mean = c(ncp_ORR, ncp_OS),
                          corr = cr))
  
}
min(power_OS)
max(loss)

