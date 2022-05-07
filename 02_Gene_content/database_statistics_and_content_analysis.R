library(svglite)
library(tidyverse)
library(data.table)
PATH <- getwd()

org.df <- fread(file.path(dirname(PATH),'01_Database_construction', 'DB_csvs', 'final_filtered.csv'))

COLORS <- c('#14f1eb', '#FFA500', '#ff0000', '#0000ff', '#008000', '#800080', '#b45f06', '#ffc0cb', '#8fce00')

grouped_org.df <- org.df %>%
  group_by(class) %>%
  summarise(mean_genes = mean(Genes),
            std_genes = sd(Genes),
            class_count = n(),
            phylum = first(phylum),
            mean_nontrna = mean(non_trna_Genes),
            std_nontrna = sd(non_trna_Genes),
            mean_trna = mean(Genes - non_trna_Genes),
            std_trna = sd(Genes - non_trna_Genes))%>%
  filter(class_count > 10) %>%
  filter(phylum != 'Placozoa') %>% 
  arrange(by_group = phylum)

svglite(file.path(PATH, 'figures', 'fig_1b_genes_per_ph.svg'), width = 6, height = 5)
g <- ggplot(data = grouped_org.df,
       aes(x = class,
           y = mean_genes,
           fill = phylum)) +
  geom_errorbar(aes(ymin = mean_genes - std_trna,
                  ymax = mean_genes + std_trna)) + 
  geom_col(alpha = .8) +
  geom_pointrange(aes(x = class, y = mean_nontrna,
                      ymin = mean_nontrna - std_nontrna,
                      ymax = mean_nontrna + std_nontrna,),
                  size = .3) +
  scale_x_discrete(limits = grouped_org.df$class) +
  scale_fill_manual(values = COLORS)+
  geom_hline(yintercept = 37, color = 'grey50', size = .7, alpha = .6) +
  geom_hline(yintercept = 15, color = 'black', size = .7, alpha = .6) +
  ylab('# Genes') +
  xlab('Class') +
  theme(axis.text = element_text(size = 12), axis.ticks = element_text(size = 12)) +
  coord_flip() + 
  theme_classic()
print(g)
dev.off()

