
library(tidyverse)
library(boot)
library(permute)
library(sampler)
library(ggthemes)
library(ggridges)
library(data.table)

# This code can be run in two different modes, with or without including tRNA genes
MODE <- 'notrna'
MODE <- 'nocodons'

COLORS <- c('#14f1eb', '#FFA500', '#ff0000', '#0000ff', '#008000', '#800080', '#b45f06', '#ffc0cb', '#8fce00')
PATH <- getwd()
org_df <- read.csv(file.path(PATH, paste('final_', MODE, '.csv', sep = '')))%>%
  filter(not_valid == 'False')
org_df$Gene_order <- sub('\\*', '', org_df$Gene_order) 

cutoff=20
calculate_class_AR_rate <- function(df, cutoff = 20){
  #Calculate the class-level AR rate for a given org data frame by first grouping by
  #class and counting the total organisms within that class and then grouping by both
  #class and gene order and counting the amount of different gene orders within each
  #class
  
  org_df <- data.frame(df)
  N_per_group <- org_df %>%
    group_by(class) %>%
    summarise(N = n(), .groups = 'drop_last') %>%
    arrange(desc(N)) %>%
    filter(N > cutoff)
  
  AR_rate_df <- org_df %>%
    group_by(class, Gene_order) %>%
    summarise(AR = n(), .groups = 'drop_last') %>%
    group_by(class) %>%
    summarise(AR = n(), .groups = 'drop_last') %>%
    inner_join(N_per_group, by = 'class') %>%
    mutate(AR_rate = ((AR - 1)/(N - 1)) * 100) %>%
    arrange(desc(AR_rate)) %>%
    left_join(select(org_df, c('class', 'phylum')), by = 'class') %>%
    distinct() %>%
    filter(class != '') %>%
    mutate(AR_rate = round(AR_rate, 1))
}

permute_class_AR_rate <- function(df, n = 10, cutoff = 20){
  # Permute the classes within my databases n times and re-calculate AR-rates for each permutation
  
  org_df <- data.frame(df)
  N <- nrow(org_df)
  
  for (i in seq(1, n)) {
    permuted <- shuffle(N)
    org_df.shuffled <- cbind(select(org_df, -Gene_order), Gene_order = org_df$Gene_order[permuted])
    AR_rate_permuted <- calculate_class_AR_rate(org_df.shuffled, cutoff = cutoff)
    
    if (i > 1) {
      total_AR <- rbind(total_AR, AR_rate_permuted)
   }else {
      total_AR <- AR_rate_permuted
   }
  }
  return(total_AR)
}

bootstrap_AR_rate <- function(df, n = 10, size = 100, cutoff = 20, permute = T)
  # Sample (with replacement) [size] organisms from each class and permute the class labels, this time on a distribution with equal sample size. This is done for [n] iterations. 
  {
  # Remove classes that are lower than cutoff
  org_df <- data.frame(df)
  classes <- org_df %>%
    count(class, sort = T) %>%
    filter(n > cutoff, class != '') %>%
    select(class)
  # Iterate [n] times (amount of bootstrap iterations)
  for (i in seq(1, n)) {
    count <- 0
    # Sample [size] organisms from each class random.y, with replacement and then add everything into a single data frame.
    for (classname in classes$class){
      booted.df <- org_df %>%
        filter(class == classname) %>%
        slice_sample(n = size, replace = TRUE)
      if (count >= 1){
        booted.total <- rbind(booted.total, booted.df)
      }
      else{
        booted.total <- booted.df
      }
      count = count + 1
    }
    # If permute mode is on, before calculating each class's AR_rates, shuffle the classes randomly.
    if (permute == T){
      N <- nrow(booted.total)
      permuted <- shuffle(N)
      booted.total.shuffled <- cbind(select(booted.total, -Gene_order), Gene_order = booted.total$Gene_order[permuted])
      AR_rate_permuted <- calculate_class_AR_rate(booted.total.shuffled)
    }
    else{
      AR_rate_permuted <- calculate_class_AR_rate(booted.total)
    }
    if (i > 1) {
      total_AR <- rbind(total_AR, AR_rate_permuted)
    }else {
      total_AR <- AR_rate_permuted
    }
  }
  return(total_AR)
}

add_p_values <- function(exp, obs, alpha = 0.05, n_test = 1, col = 'AR_rate', test){

  # Add p_values to the data frames based on a two sample test (either t-test or Mann-Whitney) between the expected and observed values for for each class. The alpha is automatically Bonferroni-corrected.
  alpha <- alpha/n_test
  obs_mean <- mean(obs)
  obs_sd <- sd(obs)
  exp_mean <- mean(exp)
  exp_sd <- sd(exp)
  if(test == 'auto'){
    exp.norm.pval <- ks.test(x = exp, 'pnorm', mean = exp_mean, sd = exp_sd)$p.value
    obs.norm.pval <- ks.test(x = obs, 'pnorm', mean = obs_mean, sd = obs_sd)$p.value
    if((exp.norm.pval < 0.05) | (obs.norm.pval < 0.05)){
      print('Not normal, performing Mann-Whitney')
      test <- 'manwt'
    }
    else{
      print('Normal, performing t-test')
      test <- 't'
    }
  }
  if(test == 'manwt'){
    p.value <- wilcox.test(x = obs, y = exp)$p.value
  }
  else if(test == 't'){
    p.value <- t.test(x = obs, y = exp)$p.value
  }
  else{
    print("Please choose either \"auto\", [manwt] or [t]! performing manwt test.")
    p.value <- wilcox.test(x = obs, y = exp)$p.value
  }
  print(paste('P.value = ', as.character(p.value)))
  
  asterisks <- case_when(
    p.value <= 0.002 * alpha ~ '****',
    p.value <= 0.02 * alpha ~ '***',
    p.value <= 0.2 * alpha ~ '**',
    p.value <= alpha ~ '*',
    p.value > alpha ~ 'ns'
  )
  return(asterisks)
}

make_asterisk_df <- function(df_obs, df_exp, alpha = 0.05, col = 'AR_rate', test = 'auto'){
  # Append asterisks to the dfs based on significance levels, for adding to plot
  asterisks <- vector()
  n_test <- length(unique(df_obs$class))
  for(i in unique((df_obs$class))){
    class_obs <- df_obs %>%
      filter(class == i) %>%
      select(AR_rate) %>% 
      as.matrix()
    
    class_exp <- df_exp %>%
      filter(class == i) %>%
      select(AR_rate) %>%
      as.matrix()
    
    class_ast <- add_p_values(class_exp, class_obs, alpha = alpha, col = col, n_test = n_test, test = test)
    asterisks <- append(asterisks, class_ast)
  }
  return(asterisks)
}

n = 10000

AR_rate_df <- calculate_class_AR_rate(org_df, cutoff = 20)

# The expected bootstrapped values are calculated by sampling 21 organisms from each class, shuffling all 21*21 organisms between classes, and calculating AR_rates.
expected <- bootstrap_AR_rate(org_df, size = 21, n = 10000, permute = TRUE)

# The observed bootstrapped values are calculated by sampling 21 organisms from each class and calculating AR_rates WITHOUT SHUFFLING.
observed <- bootstrap_AR_rate(org_df, size = 21, n = 10000, permute = FALSE)

expected <- expected %>% arrange(class)
observed <- observed %>% arrange(class)
observed$diff <- abs(observed$AR_rate - expected$AR_rate)
observed$ratio <- observed$AR_rate/expected$AR_rate


write.csv(observed, paste('observed_boot_', MODE, '.csv', sep = ''))
write.csv(expected, paste('expected_boot_', MODE, '.csv', sep = ''))

# RUN FROM HERE IF BOOTSTRAP HAS ALREADY BEEN PERFORMED
cols <- c('class', 'AR', 'AR_rate', 'phylum', 'N', 'diff', 'ratio')
observed <- fread(paste('observed_boot_', MODE, '.csv', sep = ''), select = cols)
expected <- fread(paste('expected_boot_', MODE, '.csv', sep = ''), select = cols)

ratio.median <- median(observed$ratio)
ratio.mad <- mad(observed$ratio)
asterisks <- make_asterisk_df(df_obs = observed, df_exp = expected)

boot_ratios <- vector()
N <- length(rownames(expected))
exp_copy <- copy(expected)
# This code iterates 10,000 times, for each run it randomly selects half of the expected data set twice, calculates the ratio between all the values, and adds the mean ratio to a vector. TODO - Add the std values of this vector to the graph
for(i in seq(10000)){
  first <- exp_copy %>% slice_sample(prop = .5)
  second <- exp_copy %>% slice_sample(prop = .5)
  boot_ratios <- append(boot_ratios, mean(first$AR_rate/mean(second$AR_rate)))
}


# Ratio between observed and expected
library(svglite)
x <- 0.25
svglite(paste('figures', 'fig_2d_AR_rate_bootstrapped_ratio_', MODE, '.svg', sep = ''), height = 6, width = 6)
g <- ggplot(data = observed) + 
  geom_boxplot(aes(y = ratio, x = reorder(class, AR_rate), fill = phylum), alpha = .8) +
  scale_fill_manual(values = COLORS) +
  xlab('Class') + 
  ylab('AR_rate(obs/exp)') +
  geom_hline(yintercept = 1.0, linetype = 'dashed', alpha = .5) +
  theme_classic() +
  geom_rect(data = observed[1,1], ymin = 1 - x, ymax = 1 + x, xmin = -Inf, xmax = +Inf, fill = 'grey80', alpha = .4) +
  geom_hline(yintercept = 1.0 - x, linetype = 'dashed', alpha = .7, color = 'grey80') +
  geom_hline(yintercept = 1.0 + x, linetype = 'dashed', alpha = .7, color = 'grey80') +
  coord_flip()
print(g)
dev.off()
