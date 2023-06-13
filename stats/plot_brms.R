library(ggplot2)
library(data.table)


# Full model
data <- fread("C:/git/STOC/stats/models/output_brms_0506_1926210.csv")

gg <- ggplot(data, aes(x = Trait, y = Estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = `l-95%_CI`, ymax = `u-95%_CI`))
gg
theme_set(theme_bw())

for(x in unique(data$Var)){
  df <- subset(data, data$Var == x)
  gg <- ggplot(df, aes(x = Trait, y = Estimate)) +        # ggplot2 plot with confidence intervals
    geom_point() +
    geom_pointrange(aes(ymin = `l-95%_CI`, ymax = `u-95%_CI`)) +
    ggtitle(label = x) + geom_hline(yintercept=0, linetype="dashed", color = "red") #+
    # theme(axis.text.x = element_text(angle=90))
  print(gg)
}


### Without squared term interactions  
data <- fread("C:/git/STOC/stats/models/output_brms_no2_1926212.csv")

gg <- ggplot(data, aes(x = Trait, y = Estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = `l-95%_CI`, ymax = `u-95%_CI`))
gg
theme_set(theme_bw())

for(x in unique(data$Var)){
  df <- subset(data, data$Var == x)
  gg <- ggplot(df, aes(x = Trait, y = Estimate)) +        # ggplot2 plot with confidence intervals
    geom_point() +
    geom_pointrange(aes(ymin = `l-95%_CI`, ymax = `u-95%_CI`)) +
    ggtitle(label = x) + geom_hline(yintercept=0, linetype="dashed", color = "red") #+
  # theme(axis.text.x = element_text(angle=90))
  print(gg)
}


### Only early
data <- fread("C:/git/STOC/stats/models/output_brms_early_1926218.csv")

theme_set(theme_bw())

for(x in unique(data$Var)){
  df <- subset(data, data$Var == x)
  gg <- ggplot(df, aes(x = Trait, y = Estimate)) +        # ggplot2 plot with confidence intervals
    geom_point() +
    geom_pointrange(aes(ymin = `l-95%_CI`, ymax = `u-95%_CI`)) +
    ggtitle(label = x) + geom_hline(yintercept=0, linetype="dashed", color = "red") #+
  # theme(axis.text.x = element_text(angle=90))
  print(gg)
}
