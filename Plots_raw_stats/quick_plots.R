library(ggplot2)
library(data.table)
library(cowplot)

setwd("~/repos/popgen.analysis.pipeline/Plots_raw_stats/")


# Looking at variable importance

vip = fread(input = "vip_df.csv")

gsub(x = vip$Variable, pattern = "_\\d+$", replacement = "")

vip[, Stat := gsub(x = Variable, pattern = "_\\d+$", replacement = "")]
vip[, Window := as.numeric(gsub(x = Variable, pattern = ".*_(\\d+$)", replacement = "\\1"))]


ggplot(vip, aes(factor(Window), Importance, col=model)) + geom_point() + geom_line(aes(Window, Importance, col=model)) + facet_wrap("Stat", ncol=1, scales = "free_y")



# Looking at raw data

raw = fread(input = "../data/bt_cpop.csv")

raw[, t1 := as.numeric(gsub(demography, pattern = "t1:(\\d+)_.*", replacement = "\\1"))]
raw[, s := as.numeric(gsub(demography, pattern = "t1:\\d+_s:(.+)_t2.*", replacement = "\\1"))]
raw[, t2 := as.numeric(gsub(demography, pattern = ".*_t2:(.*)$", replacement = "\\1"))]
raw[, Duration := t2]
raw[, Onset := t1-t2]
raw[, Strength := s]

# Melting the data

raw_m = melt(raw, measure.vars = 8:117, variable.name = "Statistic_Window", value.name = "Value")
raw_m[, Window := as.numeric(gsub(Statistic_Window, pattern = ".+_(\\d+)$", replacement = "\\1"))]
raw_m[, Statistic := gsub(Statistic_Window, pattern = "^(.+)_\\d+$", replacement = "\\1")]



# Points and smooth
# Const size
p1 = ggplot(raw_m[severity == 0], aes(Window, Value, col=factor(s_coef))) + geom_point(alpha=0.1) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

# Difficult bottlenecks:
p2 = ggplot(raw_m[demography == "t1:1680_s:0.05_t2:1600"], aes(Window, Value, col=factor(s_coef))) + geom_point(alpha=0.1) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"))


# The same plots without the points
# Const size
p1 = ggplot(raw_m[severity == 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

# Difficult bottleneck:
p2 = ggplot(raw_m[demography == "t1:1680_s:0.05_t2:1600"], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"))


# Plot with violins:
# Const size
p1 = ggplot(raw_m[severity == 0], aes(factor(Window), Value, fill=factor(s_coef))) + geom_violin(scale = "width") + facet_wrap("Statistic", scales = "free_y", nrow=2)

# Difficult bottleneck:
p2 = ggplot(raw_m[demography == "t1:1680_s:0.05_t2:1600"], aes(factor(Window), Value, fill=factor(s_coef))) + geom_violin(scale = "width") + facet_wrap("Statistic", scales = "free_y", nrow=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"), nrow = 2)





# Look at stats for all bottlenecks, central window 6, violin plots

# Just TajD
ggplot(raw_m[Statistic == "D" & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2)

# Just Fay&Wu's H
ggplot(raw_m[Statistic == "H" & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2)

# All stats

listOfStats = raw_m[,unique(Statistic)]

p_list = lapply(listOfStats, function (x) {
  ggplot(raw_m[Statistic == x & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2) + ggtitle(label=x)
})

plot_grid(plotlist = p_list, nrow = 2)


# Look at stats for all bottlenecks, all windows, smooth interpolation

p_list = lapply(listOfStats, function (x) {
  ggplot(raw_m[Statistic == x & severity != 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_grid(Onset~Duration + Strength) + ggtitle(label=x)
})

plot_grid(plotlist = p_list, nrow = 2)


