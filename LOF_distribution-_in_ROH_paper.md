### **LOF vs tolerated variation**

#### 1. Import LOF genotypes per individual

``` r
#Import population file
population_file <- read.table("79Chinese.id2pop.2.txt", header=FALSE, sep="\t",stringsAsFactors = TRUE)
colnames(population_file) <- c("ind","pop")
population_file$pop <- factor(population_file$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
population_file <- population_file %>% arrange(pop)
population_file <- population_file %>%
  mutate(ind = factor(ind, ind))

#Import genotypes
STOP_genotypes <- read.table("STOP_genotypes_5pop2.txt", header=TRUE, sep=c("\t"),stringsAsFactors = TRUE)
kable(STOP_genotypes[1:10,])
```

| CHROM |      POS | C_rhe_6 | C_rhe_7 | C_rhe_8 | C_rhe_9 | C_rhe_10 | C_rhe_11 | C_rhe_12 | C_rhe_13 | C_rhe_14 | C_rhe_15 | C_rhe_16 | C_rhe_45 | C_rhe_27 | C_rhe_28 | C_rhe_29 | C_rhe_30 | C_rhe_31 | C_rhe_32 | C_rhe_33 | C_rhe_34 | C_rhe_35 | C_rhe_36 | C_rhe_37 | C_rhe_38 | C_rhe_17 | C_rhe_18 | C_rhe_19 | C_rhe_20 | C_rhe_21 | C_rhe_22 | C_rhe_23 | C_rhe_24 | C_rhe_25 | C_rhe_26 | C_rhe_39 | C_rhe_40 | C_rhe_41 | C_rhe_42 | C_rhe_43 | C_rhe_44 | C_rhe_46 | C_rhe_47 | C_rhe_48 | C_rhe_49 | C_rhe_50 | C_rhe_51 | C_rhe_52 | C_rhe_53 | C_rhe_54 | C_rhe_55 | C_rhe_56 | C_rhe_57 | C_rhe_58 | C_rhe_59 | C_rhe_60 | C_rhe_61 | C_rhe_62 | C_rhe_63 | C_rhe_64 | C_rhe_65 | C_rhe_66 | C_rhe_67 | C_rhe_68 | C_rhe_69 | C_rhe_2 | C_rhe_3 | C_rhe_4 | C_rhe_5 | C_rhe_1 | C_rhe_70 | C_rhe_71 | C_rhe_74 | C_rhe_75 | C_rhe_79 | C_rhe_78 | C_rhe_76 | C_rhe_77 | C_rhe_73 | C_rhe_72 |
|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|:|
| chr2  | 27745889 | 0/0     | 0/0     | 0/0     | 0/0     | ./.      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/1      | 0/1      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 1/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | ./.     | ./.     | ./.     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 50713433 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.     | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 50713636 | 0/0     | 0/0     | 1/1     | 0/1     | 0/1      | 1/1      | 0/0      | 0/1      | 0/0      | 0/1      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 1/1      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/1      | 0/0      | 0/1      | 0/1      | 0/1      | 0/0      | 0/0      | 0/1      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/1      | 0/1     | 0/0     | 0/0     | 0/0     | 0/1     | 0/0      | 0/0      | 0/1      | 0/1      | 0/0      | 0/0      | 0/0      | 0/1      | 0/1      | 0/1      |
| chr2  | 50713683 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0     | 0/0     | 0/0     | 0/0     | 0/0     | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      |
| chr2  | 50713927 | 0/1     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/1      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0     | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      |
| chr2  | 57956388 | 0/0     | 0/0     | 0/1     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0     | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | ./.      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 84368555 | 0/0     | ./.     | 0/0     | 0/0     | 1/1      | ./.      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 1/1      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 1/1      | 0/1      | 0/1      | 0/1      | 0/1      | 0/1      | 0/0      | 0/0      | ./.      | 0/0      | 0/1      | 0/0      | 0/1      | 0/0      | 1/1      | 0/1      | 0/1      | 0/0      | 0/1      | 0/0      | 0/1      | 1/1      | 0/0      | 0/1      | 0/1      | 1/1      | 0/1      | 0/0      | 0/1      | 0/1      | 0/1      | 0/0      | 1/1      | 0/0      | 0/1      | 1/1      | 0/1      | 0/1      | 0/0     | ./.     | 0/0     | ./.     | 0/1     | 0/1      | 0/1      | 0/0      | 0/1      | 0/0      | 0/1      | 1/1      | 0/0      | 0/1      | 0/1      |
| chr3  |   541039 | 0/0     | 0/0     | ./.     | 0/0     | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | ./.      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0     | ./.     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |
| chr3  | 33007093 | ./.     | 0/0     | ./.     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | ./.     | 0/0     | ./.     | ./.     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |
| chr3  | 40003876 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/1      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0     | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      | 0/0      |

#### 2. Counting the number of genome-wide LOF per individual

``` r
#counting the genotypes
genotypes_count_STOP <- setNames(cbind.data.frame(rep("stop",79), as.data.frame(colnames(STOP_genotypes)[3:81])),c("pred","ind"))
genotypes_count_STOP$homoref <- colSums(STOP_genotypes[,3:ncol(STOP_genotypes)] == "0/0", na.rm = TRUE)
genotypes_count_STOP$hetero <- colSums(STOP_genotypes[,3:ncol(STOP_genotypes)] == "0/1", na.rm = TRUE)
genotypes_count_STOP$homoalt <- colSums(STOP_genotypes[,3:ncol(STOP_genotypes)] == "1/1", na.rm = TRUE)
genotypes_count_STOP$miss <- colSums(STOP_genotypes[,3:ncol(STOP_genotypes)] == "./.", na.rm = TRUE)

#adding the population name
genotypes_count_STOP$pop <- NA
for (i in 1:nrow(genotypes_count_STOP)) {
  ind=as.character(genotypes_count_STOP$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  genotypes_count_STOP[genotypes_count_STOP$ind==ind, "pop"] <- pop
}

genotypes_count_STOP <- as.data.frame(unclass(genotypes_count_STOP), stringsAsFactors = TRUE)
colnames(genotypes_count_STOP) <- c("pred","ind","homoref","hetero","homoalt","miss","pop")
genotypes_count_STOP <- genotypes_count_STOP %>%
    mutate(total = select(., homoref:miss) %>% rowSums(na.rm = TRUE))

genotypes_count_STOP$pop <- factor(genotypes_count_STOP$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
genotypes_count_STOP <- genotypes_count_STOP %>% arrange(pop)
genotypes_count_STOP <- genotypes_count_STOP %>%
  mutate(pop, ind = factor(ind, ind))

#print
kable(genotypes_count_STOP[1:10,], col.names = c("Prediction","Individual", "0/0","0/1","1/1","./.","Population", "Total"))
```

| Prediction | Individual | 0/0 | 0/1 | 1/1 | ./. | Population | Total |
|:-----------|:-----------|----:|----:|----:|----:|:-----------|------:|
| stop       | C_rhe_70   | 129 |  26 |   4 |  20 | mulatta    |   179 |
| stop       | C_rhe_71   | 144 |  23 |   7 |   5 | mulatta    |   179 |
| stop       | C_rhe_74   | 146 |  28 |   4 |   1 | mulatta    |   179 |
| stop       | C_rhe_75   | 140 |  27 |   6 |   6 | mulatta    |   179 |
| stop       | C_rhe_79   | 145 |  26 |   6 |   2 | mulatta    |   179 |
| stop       | C_rhe_78   | 133 |  32 |   5 |   9 | mulatta    |   179 |
| stop       | C_rhe_76   | 140 |  27 |   4 |   8 | mulatta    |   179 |
| stop       | C_rhe_77   | 142 |  24 |   5 |   8 | mulatta    |   179 |
| stop       | C_rhe_73   | 152 |  20 |   6 |   1 | mulatta    |   179 |
| stop       | C_rhe_72   | 145 |  27 |   4 |   3 | mulatta    |   179 |

*Figure 5*

``` r
leg <- c("mulatta" = "chartreuse3","lasiotis" = "orchid3", "brevicaudus" = "darkorange", "littoralis" = "red2", "tcheliensis" = "blue3")

#LOF HETEROZYGOTES
lofhetero <- ggplot(genotypes_count_STOP) +
  geom_violin(aes(x=pop, y=(hetero), fill = pop, color=pop), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(aes(x=pop, y=(hetero), fill = pop, color=pop), show.legend = FALSE) +
  labs(x = "Individuals", y = "Number of LOF heterozygotes", title = "B") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  guides(fill=guide_legend(title="Population")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
#LOF HOMOZYGOTES
lofhomo <- ggplot(genotypes_count_STOP) +
  geom_violin(aes(x=pop, y=(homoalt), fill = pop, color=pop), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(aes(x=pop, y=(homoalt), fill = pop, color=pop), show.legend = FALSE) +
  labs(x = "Individuals", y = "Number of LOF homozygotes", title = "A") +
  scale_colour_manual(values = leg) +
  scale_y_continuous(breaks = seq(0, 12, 2), limits = c(0, 12)) +
  scale_fill_manual(values = alpha(leg, .3)) +
  guides(fill=guide_legend(title="Population")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

#LOF ALLELES
lofallele <- ggplot(genotypes_count_STOP) +
  geom_violin(aes(x=pop, y=(2*homoalt+hetero), fill = pop, color=pop), draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) +
  geom_point(aes(x=pop, y=(2*homoalt+hetero), fill = pop, color=pop), show.legend = FALSE) +
  labs(x = "Individuals", y = "Number of LOF alleles", title = "C") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  guides(fill=guide_legend(title="Population")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

TOTLOF <- grid.arrange(lofhomo, lofhetero, lofallele, ncol=3)
```

![](LOF_distribution-_in_ROH_paper_files/figure-markdown_github/stophist-1.png)

``` r
#ggsave("Total_LOF_variation.png", plot = TOTLOF, width = 30, height = 11, units = "cm")
#cairo_ps(file = "Total_LOF_variation.eps", onefile = FALSE, width=11.81, height=4.33)
#plot(TOTLOF)
#dev.off()
```

Wilcox-test, two-sided

    ##       [,1]          [,2]          [,3]               
    ##  [1,] "mulatta"     "mulatta"     "1"                
    ##  [2,] "mulatta"     "lasiotis"    "0.444726613344227"
    ##  [3,] "mulatta"     "littoralis"  "0.76264881499372" 
    ##  [4,] "mulatta"     "brevicaudus" "0.801206422946389"
    ##  [5,] "mulatta"     "tcheliensis" "0.23650704969481" 
    ##  [6,] "lasiotis"    "mulatta"     "0.444726613344227"
    ##  [7,] "lasiotis"    "lasiotis"    "1"                
    ##  [8,] "lasiotis"    "littoralis"  "0.262367731839683"
    ##  [9,] "lasiotis"    "brevicaudus" "0.833431842167942"
    ## [10,] "lasiotis"    "tcheliensis" "0.1428282927079"  
    ## [11,] "littoralis"  "mulatta"     "0.76264881499372" 
    ## [12,] "littoralis"  "lasiotis"    "0.262367731839683"
    ## [13,] "littoralis"  "littoralis"  "1"                
    ## [14,] "littoralis"  "brevicaudus" "0.593417622733999"
    ## [15,] "littoralis"  "tcheliensis" "0.59367061762224" 
    ## [16,] "brevicaudus" "mulatta"     "0.801206422946389"
    ## [17,] "brevicaudus" "lasiotis"    "0.833431842167942"
    ## [18,] "brevicaudus" "littoralis"  "0.593417622733999"
    ## [19,] "brevicaudus" "brevicaudus" "1"                
    ## [20,] "brevicaudus" "tcheliensis" "0.39166779226774" 
    ## [21,] "tcheliensis" "mulatta"     "0.23650704969481" 
    ## [22,] "tcheliensis" "lasiotis"    "0.1428282927079"  
    ## [23,] "tcheliensis" "littoralis"  "0.59367061762224" 
    ## [24,] "tcheliensis" "brevicaudus" "0.39166779226774" 
    ## [25,] "tcheliensis" "tcheliensis" "1"

Wilcox-test, less

    ##       [,1]          [,2]          [,3]                 
    ##  [1,] "mulatta"     "mulatta"     "0.515366481925404"  
    ##  [2,] "mulatta"     "lasiotis"    "0.642739370326541"  
    ##  [3,] "mulatta"     "littoralis"  "0.986658488448782"  
    ##  [4,] "mulatta"     "brevicaudus" "0.976099091458559"  
    ##  [5,] "mulatta"     "tcheliensis" "0.990401850948536"  
    ##  [6,] "lasiotis"    "mulatta"     "0.368696326055292"  
    ##  [7,] "lasiotis"    "lasiotis"    "0.502825337161205"  
    ##  [8,] "lasiotis"    "littoralis"  "0.994680477779723"  
    ##  [9,] "lasiotis"    "brevicaudus" "0.932522408168885"  
    ## [10,] "lasiotis"    "tcheliensis" "0.987151718695554"  
    ## [11,] "littoralis"  "mulatta"     "0.0145253204699465" 
    ## [12,] "littoralis"  "lasiotis"    "0.00555711655430629"
    ## [13,] "littoralis"  "littoralis"  "0.503293292746948"  
    ## [14,] "littoralis"  "brevicaudus" "0.52015741850953"   
    ## [15,] "littoralis"  "tcheliensis" "0.867250166120553"  
    ## [16,] "brevicaudus" "mulatta"     "0.0317673759265758" 
    ## [17,] "brevicaudus" "lasiotis"    "0.0736899183117097" 
    ## [18,] "brevicaudus" "littoralis"  "0.5"                
    ## [19,] "brevicaudus" "brevicaudus" "0.542235013311614"  
    ## [20,] "brevicaudus" "tcheliensis" "0.827866725954459"  
    ## [21,] "tcheliensis" "mulatta"     "0.0132620934481365" 
    ## [22,] "tcheliensis" "lasiotis"    "0.014452788850682"  
    ## [23,] "tcheliensis" "littoralis"  "0.143919729012942"  
    ## [24,] "tcheliensis" "brevicaudus" "0.230987715519692"  
    ## [25,] "tcheliensis" "tcheliensis" "0.543327939942697"

Wilcox-test, greater

    ##       [,1]          [,2]          [,3]                 
    ##  [1,] "mulatta"     "mulatta"     "0.515224485176494"  
    ##  [2,] "mulatta"     "lasiotis"    "0.65436369686717"   
    ##  [3,] "mulatta"     "brevicaudus" "0.229603658377577"  
    ##  [4,] "mulatta"     "littoralis"  "0.0363504832167738" 
    ##  [5,] "mulatta"     "tcheliensis" "0.0785655112564528" 
    ##  [6,] "lasiotis"    "mulatta"     "0.356966061235359"  
    ##  [7,] "lasiotis"    "lasiotis"    "0.502826490105507"  
    ##  [8,] "lasiotis"    "brevicaudus" "0.134946149087358"  
    ##  [9,] "lasiotis"    "littoralis"  "0.00272574401912657"
    ## [10,] "lasiotis"    "tcheliensis" "0.0313045275559253" 
    ## [11,] "brevicaudus" "mulatta"     "0.806071860401777"  
    ## [12,] "brevicaudus" "lasiotis"    "0.874780118657248"  
    ## [13,] "brevicaudus" "brevicaudus" "0.542235013311614"  
    ## [14,] "brevicaudus" "littoralis"  "0.334338790189806"  
    ## [15,] "brevicaudus" "tcheliensis" "0.273809523809524"  
    ## [16,] "littoralis"  "mulatta"     "0.966220621080591"  
    ## [17,] "littoralis"  "lasiotis"    "0.997399367002806"  
    ## [18,] "littoralis"  "brevicaudus" "0.683786036316919"  
    ## [19,] "littoralis"  "littoralis"  "0.503281022381286"  
    ## [20,] "littoralis"  "tcheliensis" "0.157095714788737"  
    ## [21,] "tcheliensis" "mulatta"     "0.937952670585874"  
    ## [22,] "tcheliensis" "lasiotis"    "0.971799993782548"  
    ## [23,] "tcheliensis" "brevicaudus" "0.78968253968254"   
    ## [24,] "tcheliensis" "littoralis"  "0.854696135294011"  
    ## [25,] "tcheliensis" "tcheliensis" "0.542235013311614"

#### 3. LOF vs tolerated homozygotes in ROH

``` r
stop_summary <- read.table("STOP_distribution_5pop.txt", header = TRUE, stringsAsFactors = TRUE)
tol_summary <- read.table("tolerated_distribution_5pop.txt", header = TRUE, stringsAsFactors = TRUE)
summary_tolstop <- rbind(setNames(cbind(rep("stop", 79), stop_summary),c("pred",colnames(stop_summary))), setNames(cbind(rep("tolerated", 79), tol_summary),c("pred",colnames(stop_summary))))

options(digits=3)
pop <- c("mulatta","lasiotis", "brevicaudus", "littoralis", "tcheliensis")

summary_tolstop[["deleterious"]] <- as.numeric(summary_tolstop$pred=='stop')
summary_tolstop[["frac.inA"]] <- summary_tolstop$nA/summary_tolstop$ntotal
summary_tolstop[["frac.inB"]] <- summary_tolstop$nB/summary_tolstop$ntotal
summary_tolstop[["frac.inC"]] <- summary_tolstop$nC/summary_tolstop$ntotal
summary_tolstop[["frac.inany"]] <- summary_tolstop$nallroh/summary_tolstop$ntotal

#adding the population name
summary_tolstop$pop <- NA
for (i in 1:nrow(summary_tolstop)) {
  ind=as.character(summary_tolstop$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  summary_tolstop[summary_tolstop$ind==ind, "pop"] <- pop
}

#write.csv(summary_tolstop, file = "stoptol_distributionROH.csv")
kable(summary_tolstop[1:10,])
```

| pred | ind      |  nA |  nB |  nC | nallroh | ntotal | deleterious | frac.inA | frac.inB | frac.inC | frac.inany | pop         |
|:---|:------|--:|--:|--:|-----:|-----:|-------:|------:|------:|------:|-------:|:-------|
| stop | C_rhe_1  |   0 |   2 |   2 |       4 |      9 |           1 |    0.000 |    0.222 |    0.222 |      0.444 | tcheliensis |
| stop | C_rhe_10 |   0 |   0 |   1 |       1 |      5 |           1 |    0.000 |    0.000 |    0.200 |      0.200 | littoralis  |
| stop | C_rhe_11 |   0 |   0 |   0 |       0 |      7 |           1 |    0.000 |    0.000 |    0.000 |      0.000 | littoralis  |
| stop | C_rhe_12 |   0 |   0 |   0 |       0 |      0 |           1 |      NaN |      NaN |      NaN |        NaN | littoralis  |
| stop | C_rhe_13 |   0 |   0 |   1 |       1 |      7 |           1 |    0.000 |    0.000 |    0.143 |      0.143 | littoralis  |
| stop | C_rhe_14 |   0 |   0 |   1 |       1 |      4 |           1 |    0.000 |    0.000 |    0.250 |      0.250 | littoralis  |
| stop | C_rhe_15 |   0 |   0 |   1 |       1 |      7 |           1 |    0.000 |    0.000 |    0.143 |      0.143 | littoralis  |
| stop | C_rhe_16 |   0 |   0 |   0 |       0 |      5 |           1 |    0.000 |    0.000 |    0.000 |      0.000 | littoralis  |
| stop | C_rhe_17 |   1 |   1 |   1 |       3 |      8 |           1 |    0.125 |    0.125 |    0.125 |      0.375 | littoralis  |
| stop | C_rhe_18 |   1 |   0 |   0 |       1 |      3 |           1 |    0.333 |    0.000 |    0.000 |      0.333 | littoralis  |

#### 4. Importing ROH coverage per individual

``` r
roh <- read.table("5pop_roh_summary.txt", header=TRUE, sep="\ ", stringsAsFactors = TRUE)
roh[["ind"]] <- rownames(roh)
rownames(roh) <- NULL
summary_tolstop[["cov.any"]] <- NA 
summary_tolstop[["cov.A"]] <- NA
summary_tolstop[["cov.B"]] <- NA
summary_tolstop[["cov.C"]] <- NA

genome_size=2900000000

for (i in 1:nrow(summary_tolstop)) {
  ind=summary_tolstop[i,]$ind
  summary_tolstop[i,"cov.A"] <- roh[roh$ind==ind,][1,"A"]/genome_size
  summary_tolstop[i,"cov.B"] <- roh[roh$ind==ind,][1,"B"]/genome_size
  summary_tolstop[i,"cov.C"] <- roh[roh$ind==ind,][1,"C"]/genome_size
  summary_tolstop[i,"cov.any"] <- roh[roh$ind==ind,][1,"TOTAL"]/genome_size
}

summary_tolstop[["highcovA"]] <- as.numeric(summary_tolstop$cov.A>=quantile(summary_tolstop$cov.A)[3])
summary_tolstop[["highcovB"]] <- as.numeric(summary_tolstop$cov.B>=quantile(summary_tolstop$cov.B)[3])
summary_tolstop[["highcovC"]] <- as.numeric(summary_tolstop$cov.C>=quantile(summary_tolstop$cov.C)[3])
summary_tolstop[["highcovany"]] <- as.numeric(summary_tolstop$cov.any>=quantile(summary_tolstop$cov.any)[3])

summary_tolstop[,c(18,19,20,21)][summary_tolstop[,c(18,19,20,21)]==1] <- "high"
summary_tolstop[,c(18,19,20,21)][summary_tolstop[,c(18,19,20,21)]==0] <- "low"

resume_tolstop <- summary_tolstop[,c(1,2,9,10,11,12,13,18,19,20,21)]
resume_tolstop <- as.data.frame(unclass(resume_tolstop), stringsAsFactors = TRUE)
resume_tolstop$pred <- factor(resume_tolstop$pred, level=c("tolerated", "stop"))
resume_tolstop$highcovA <- factor(resume_tolstop$highcovA, level=c("low", "high"))
resume_tolstop$highcovB <- factor(resume_tolstop$highcovB, level=c("low", "high"))
resume_tolstop$highcovC <- factor(resume_tolstop$highcovC, level=c("low", "high"))
resume_tolstop$highcovany <- factor(resume_tolstop$highcovany, level=c("low", "high"))

options(digits=3)
#write.csv(summary_tolstop, file = "data_stop_tol_inROH.csv", row.names = FALSE, quote=FALSE)
kable(resume_tolstop[1:10,])
```

| pred | ind      | frac.inA | frac.inB | frac.inC | frac.inany | pop         | highcovA | highcovB | highcovC | highcovany |
|:---|:------|------:|------:|------:|-------:|:--------|:------|:------|:------|:-------|
| stop | C_rhe_1  |    0.000 |    0.222 |    0.222 |      0.444 | tcheliensis | high     | high     | high     | high       |
| stop | C_rhe_10 |    0.000 |    0.000 |    0.200 |      0.200 | littoralis  | high     | high     | low      | high       |
| stop | C_rhe_11 |    0.000 |    0.000 |    0.000 |      0.000 | littoralis  | high     | high     | high     | high       |
| stop | C_rhe_12 |      NaN |      NaN |      NaN |        NaN | littoralis  | high     | high     | high     | high       |
| stop | C_rhe_13 |    0.000 |    0.000 |    0.143 |      0.143 | littoralis  | high     | high     | high     | high       |
| stop | C_rhe_14 |    0.000 |    0.000 |    0.250 |      0.250 | littoralis  | high     | high     | high     | high       |
| stop | C_rhe_15 |    0.000 |    0.000 |    0.143 |      0.143 | littoralis  | high     | high     | low      | low        |
| stop | C_rhe_16 |    0.000 |    0.000 |    0.000 |      0.000 | littoralis  | high     | high     | high     | high       |
| stop | C_rhe_17 |    0.125 |    0.125 |    0.125 |      0.375 | littoralis  | low      | low      | low      | low        |
| stop | C_rhe_18 |    0.333 |    0.000 |    0.000 |      0.333 | littoralis  | low      | low      | low      | low        |

``` r
##TESTS
options(digits = 3)

#Wilcoxon tests TWO-SIDED
wilcoxtest_anyhigh2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inany, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inany)
wilcoxtest_anylow2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inany, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inany)
wilcoxtest_Ahigh2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inA, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inA)
wilcoxtest_Alow2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inA, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inA)
wilcoxtest_Bhigh2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inB, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inB)
wilcoxtest_Blow2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inB, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inB)
wilcoxtest_Chigh2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inC, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inC)
wilcoxtest_Clow2s <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inC, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inC)

#Wilcoxon tests ONE SIDED LOF < TOL
wilcoxtest_anyhigh <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inany, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inany, alternative = "less")
wilcoxtest_anylow <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inany, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inany, alternative = "less")
wilcoxtest_Ahigh <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inA, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inA, alternative = "less")
wilcoxtest_Alow <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inA, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inA, alternative = "less")
wilcoxtest_Bhigh <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inB, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inB, alternative = "less")
wilcoxtest_Blow <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inB, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inB, alternative = "less")
wilcoxtest_Chigh <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="stop", ]$frac.inC, resume_tolstop[resume_tolstop$highcovany=="high" & resume_tolstop$pred=="tolerated", ]$frac.inC, alternative = "less")
wilcoxtest_Clow <- wilcox.test(resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="stop", ]$frac.inC, resume_tolstop[resume_tolstop$highcovany=="low" & resume_tolstop$pred=="tolerated", ]$frac.inC, alternative = "less")
```

**Wilcox test - two-sided**

| ROH class | ROH coverage | pvalue               | W   |
|:----------|:-------------|:---------------------|:----|
| AnyROH    | Low          | 0.0122333344330271   | 513 |
| AnyROH    | High         | 0.210585965083396    | 652 |
| A         | Low          | 4.56007790734886e-07 | 273 |
| A         | High         | 3.92565701097994e-08 | 240 |
| B         | Low          | 4.56007790734886e-07 | 273 |
| B         | High         | 3.03195714086971e-06 | 318 |
| C         | Low          | 9.05755814251273e-09 | 222 |
| C         | High         | 0.41048623163095     | 696 |

**Wilcox test - LOF less than tolerated**

| ROH class | ROH coverage | pvalue               | W   |
|:----------|:-------------|:---------------------|:----|
| AnyROH    | Low          | 0.00611666721651356  | 513 |
| AnyROH    | High         | 0.105292982541698    | 652 |
| A         | Low          | 2.28003895367443e-07 | 273 |
| A         | High         | 1.96282850548997e-08 | 240 |
| B         | Low          | 2.28003895367443e-07 | 273 |
| B         | High         | 1.51597857043486e-06 | 318 |
| C         | Low          | 4.52877907125637e-09 | 222 |
| C         | High         | 0.205243115815475    | 696 |

*Figure 6*

``` r
leg <- c("stop" = "red","tolerated" = "blue")

plot_list_stop <- list()

label.df <- data.frame(highcovany = c("low"), frac.inany = c(0.8))

plot_list_stop[[1]] <- ggplot_gtable(ggplot_build(ggplot() +
  geom_boxplot(data = resume_tolstop, mapping = aes(x=highcovany, y=frac.inany, fill=pred, color=pred), outlier.shape = NA) +
  geom_point(data = resume_tolstop, mapping = aes(x=highcovany, y=frac.inany, fill=pred, color=pred), position=position_dodge(width=0.75), size = 0.65) +
  geom_text(data = label.df, mapping = aes(x=highcovany, y=frac.inany), label = c(paste("p=", signif(wilcoxtest_anylow$p.value, digits = 3), sep = "")), size = 3) +
  labs(x = "ROH density", y = paste(expression("Fraction of homozygotes\nfalling into any size ROH")), title="A") +
  scale_fill_manual(values=alpha(leg, .3), labels = c("Tolerated", "LOF")) +
  scale_color_manual(values=leg, labels = c("Tolerated","LOF")) +
  guides(fill=guide_legend(nrow = 1), color = guide_legend(nrow = 1)) +
  scale_x_discrete(labels=c("Low ROH", "High ROH")) +
  scale_y_continuous(limits = c(0, 1)) +
  #scale_colour_manual(values = leg, name="SIFT Prediction", labels = c("Tolerated", "LOF")) +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=10,color = "black"), axis.title.y=element_text(size=11), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=11), legend.position = c(0.32, 0.93), legend.text = element_text(size=9, margin = unit(c(0,0,0,0.15), "cm")), legend.key.height = unit(0.3, 'cm'), legend.title =  element_blank(), legend.key.spacing.x = unit(0.5, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), axis.ticks = element_line(colour = "black")))) 


label.df <- data.frame(highcovA = c("high", "low"), frac.inA = c(0.4, 0.4))


plot_list_stop[[2]] <- ggplot_gtable(ggplot_build(ggplot() +
  geom_boxplot(data = resume_tolstop, mapping = aes(x=highcovA, y=frac.inA, fill=pred, color=pred), show.legend = FALSE, outlier.shape = NA) +
  geom_point(data = resume_tolstop, mapping = aes(x=highcovA, y=frac.inA, fill=pred, color=pred), show.legend = FALSE, position=position_dodge(width=0.75), size = 0.65) +
  geom_text(data = label.df, mapping = aes(x=highcovA, y=frac.inA), label = c(paste("p=", signif(wilcoxtest_Ahigh$p.value, digits = 3), sep = ""), paste("p=", signif(wilcoxtest_Alow$p.value, digits = 3), sep = "")), size = 3) +
  labs(x = "ROH density", y = paste(expression("Fraction of homozygotes\nfalling in short ROH", sep="")), title="B") +
  scale_fill_manual(values=alpha(leg, .35)) +
  scale_color_manual(values=leg) +
  guides(fill=guide_legend(title="SIFT Prediction"), colour=guide_legend(title="SIFT Prediction")) +
  scale_x_discrete(labels=c("Low ROH", "High ROH")) +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=10,color = "black"), axis.title.y=element_text(size=11), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=13), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), axis.ticks = element_line(colour = "black"))))

label.df <- data.frame(highcovB = c("high", "low"), frac.inB = c(0.3, 0.3))

plot_list_stop[[3]] <- ggplot_gtable(ggplot_build(ggplot() +
  geom_boxplot(data = resume_tolstop, mapping = aes(x=highcovB, y=frac.inB, fill=pred, color=pred), show.legend = FALSE, outlier.shape = NA) +
  geom_point(data = resume_tolstop, mapping = aes(x=highcovB, y=frac.inB, fill=pred, color=pred), show.legend = FALSE, position=position_dodge(width=0.75), size = 0.65) +
  geom_text(data = label.df, mapping = aes(x=highcovB, y=frac.inB), label = c(paste("p=", formatC(wilcoxtest_Bhigh$p.value, format = "e", digits = 2), sep = ""), paste("p=", signif(wilcoxtest_Blow$p.value, digits = 3), sep = "")), size = 3) +
  labs(x = "ROH density", y = paste(expression("Fraction of homozygotes\nfalling in medium ROH", sep="")), title="C") +
  scale_fill_manual(values=alpha(leg, .3)) +
  scale_color_manual(values=leg) +
  guides(fill=guide_legend(title="SIFT Prediction"), colour=guide_legend(title="SIFT Prediction")) +
  scale_x_discrete(labels=c("Low ROH", "High ROH")) +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=10,color = "black"), axis.title.y=element_text(size=11), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=13), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), axis.ticks = element_line(colour = "black"))))

label.df <- data.frame(highcovC = c("low"), frac.inC = c(0.75))

plot_list_stop[[4]] <- ggplot_gtable(ggplot_build(ggplot() +
  geom_boxplot(data = resume_tolstop, mapping = aes(x=highcovC, y=frac.inC, fill=pred, color=pred), show.legend = FALSE, outlier.shape = NA) +
  geom_point(data = resume_tolstop, mapping = aes(x=highcovC, y=frac.inC, fill=pred, color=pred), show.legend = FALSE, position=position_dodge(width=0.75), size = 0.65) +
  geom_text(data = label.df, mapping = aes(x=highcovC, y=frac.inC), label = c(paste("p=", signif(wilcoxtest_Clow$p.value, digits = 2), sep = "")), size = 3) +
  labs(x = "ROH density", y = paste(expression("Fraction of homozygotes\nfalling in long ROH", sep="")), title="D") +
  scale_fill_manual(values=alpha(leg, .3)) +
  scale_color_manual(values=leg) +
  guides(fill=guide_legend(title="SIFT Prediction"), colour=guide_legend(title="SIFT Prediction")) +
  scale_x_discrete(labels=c("Low ROH", "High ROH")) +
  theme_light() +
  theme(axis.title.x = element_blank(), axis.text = element_text(size=10,color = "black"), axis.title.y=element_text(size=11), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=13), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), axis.ticks = element_line(colour = "black"))))

fig <- do.call("grid.arrange", c(plot_list_stop, ncol=2))
```

![](LOF_distribution-_in_ROH_paper_files/figure-markdown_github/stopplot-1.png)

``` r
#ggsave("Frac_LOF_ROH.png", fig, width = 20, height = 15, units = "cm")

#cairo_ps(file = "Frac_LOF_ROH.eps", onefile = FALSE, width=7.87, height=5.9)
#plot(fig)
#dev.off()
```
