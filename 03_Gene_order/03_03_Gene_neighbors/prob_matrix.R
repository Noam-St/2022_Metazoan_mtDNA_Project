library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(Cairo)
library(gmodels)

prev_long <- read.csv(file.path('prob_matrices', 'combined_prob_long_notrna.csv'))
ss_long <- read.csv(file.path('prob_matrices', 'combined_prob_long_notrna.csv'))
complex1 <- c('nad1', 'nad2', 'nad3', 'nad4', 'nad4L', 'nad5', 'nad6')
complex3 <- c('cob')
complex4 <- c('cox1', 'cox2', 'cox3')
complex5 <- c('atp6', 'atp8')
rrna <- c('rnl', 'rns')

prev_long %>%
  pivot_longer(!c(gene_1, gene_2), names_to = 'Phylum', values_to = 'Prevalence') %>%
  ggplot(aes(x = gene_1, y = gene_2, fill = Prevalence)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), axis.text.y = element_text(vjust = .5, hjust = 1)) + 
  xlab('') + 
  ylab('') +
  facet_wrap(vars(Phylum))

prev_long_nodups <- prev_long %>% 
  filter(gene_1 != gene_2) %>%
  rowwise() %>%
  mutate(Mean = mean(c(Arthropoda, Chordata, Cnidaria, Echinodermata, Mollusca))) %>%
  filter(Mean > 0.05) %>% 
  select(-Mean)
ss_long_mat <- ss_long %>% 
  semi_join(prev_long_nodups, by = c('gene_1', 'gene_2')) %>% 
  select(-c(gene_1, gene_2)) %>% 
  as.matrix()

heat_data <- prev_long_nodups %>%
  select(-c(gene_1, gene_2)) %>%
  as.matrix() 
rownames(heat_data) <- paste(prev_long_nodups$gene_1, prev_long_nodups$gene_2)
pheatmap(heat_data, cluster_rows = T, cutree_rows = 3, border_color = 'grey60', show_rownames = T, fontsize = 7.5, height = 6.5, width = 4, filename = 'heatmap.png')

col_fun <- colorRamp2(c(0, 0.5, 1), c('white','orange', 'red'), space = 'RGB')
custom_legend <- Legend(at = '',
                        legend_gp = gpar(fill = 'black'),
                        title = 'Different strand',
                        type = 'points')
map <- Heatmap(heat_data,
        name = 'Prevalence',
        col = col_fun,
        border_gp = gpar(col = 'black', lty = 1),
        rect_gp = gpar(col = 'grey90'),
        column_title = 'Phylum',
        column_title_side = 'bottom',
        show_row_names = T,
        row_names_gp = gpar(fontsize = 7,'black'),
        row_names_centered = F,
        column_names_centered = T,
        column_names_rot = 45,
        use_raster = T,
        row_km = 3,
        row_gap = unit(1.5, 'mm'),
        cell_fun = function(j, i, x, y, width, height, fill){
          if(ss_long_mat[i, j] == 'False' ){
            grid.circle(x = x,
                        y = y,
                        r = unit(.08, 'cm'),
                        gp = gpar(fill = 'black', col = 'black'))
          }
        })
draw(map, annotation_legend_list = custom_legend, ht_gap = unit(7, 'mm'),heatmap_legend_side = 'right', annotation_legend_side = 'right' )

ss_longer <- ss_long %>%
  semi_join(prev_long_nodups, by = c('gene_1', 'gene_2')) %>% 
  pivot_longer(!c(gene_1, gene_2), names_to = 'Phylum', values_to = 'Same_strand')
  
prev_long_nodups_longer <- prev_long_nodups %>% 
  pivot_longer(!c(gene_1, gene_2), names_to = 'Phylum', values_to = 'Prevalence')%>%
  mutate(conserved = Prevalence > .5) %>% 
  left_join(ss_longer)
sum((prev_long_nodups_longer$conserved == T) & (prev_long_nodups_longer$Same_strand == 'True'))/sum(prev_long_nodups_longer$Same_strand == 'True')
sum((prev_long_nodups_longer$conserved == T) & (prev_long_nodups_longer$Same_strand == 'False'))/sum(prev_long_nodups_longer$Same_strand == 'False')

test <- CrossTable(
  prev_long_nodups_longer$conserved,
  prev_long_nodups_longer$Same_strand,
  chisq = T,
  fisher = T,
  mcnemar = T,
  resid = T,
  asresid = T)

