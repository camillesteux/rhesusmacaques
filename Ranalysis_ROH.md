``` r
#Importing packages
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(ggpubr)
library(grid)
library(dplyr)
library(knitr)
```

### **ROH distribution analysis**

#### Importing data (ROH length and number per individuals)

``` r
#Importing data
data <- read.table("5pop_roh_summary.txt", header=TRUE, sep="\ ",stringsAsFactors = TRUE)
kable(data[1:10,])
```

|          | pop        |   nA |        A |  nB |        B |  nC |         C | nTOTAL |     TOTAL |
|:--------|:---------|----:|--------:|----:|--------:|---:|--------:|------:|--------:|
| C_rhe_6  | littoralis |  558 | 22524127 | 124 | 18533330 |   2 |    940136 |    684 |  41997593 |
| C_rhe_7  | littoralis | 1056 | 49013003 | 432 | 74766235 |  45 |  65619603 |   1533 | 189398841 |
| C_rhe_8  | littoralis | 1125 | 53157254 | 451 | 77212368 |  43 |  44225569 |   1619 | 174595191 |
| C_rhe_9  | littoralis | 1055 | 48493391 | 451 | 80009353 |  82 | 152473535 |   1588 | 280976279 |
| C_rhe_10 | littoralis | 1140 | 52702275 | 412 | 71920412 |  46 |  76559815 |   1598 | 201182502 |
| C_rhe_11 | littoralis | 1080 | 51632918 | 427 | 74644346 |  94 | 186742429 |   1601 | 313019693 |
| C_rhe_12 | littoralis |  886 | 39332213 | 283 | 50307649 |  71 | 183948530 |   1240 | 273588392 |
| C_rhe_13 | littoralis |  887 | 40367821 | 287 | 50594287 |  70 | 163745816 |   1244 | 254707924 |
| C_rhe_14 | littoralis |  898 | 39896497 | 271 | 46804620 |  88 | 205700375 |   1257 | 292401492 |
| C_rhe_15 | littoralis |  961 | 42551381 | 318 | 57541961 |  54 |  78609089 |   1333 | 178702431 |

#### Plotting data

##### 1. SROH distribution per population (Figure S1)

``` r
#Setting legend colors so it matches Liu. et al (2017) paper
leg <- c("mulatta" = "chartreuse3","lasiotis" = "orchid3", "brevicaudus" = "darkorange", "littoralis" = "red2", "tcheliensis" = "blue3")
pop_order <- factor(data$pop, level=c('mulatta', 'lasiotis', 'brevicaudus', 'littoralis','tcheliensis'))

#Plots
plotAll1 <- ggplot(data = data, mapping=aes(x=pop_order, y=TOTAL/(1e6), colour=pop_order, fill=pop_order)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) + 
  geom_point(show.legend = FALSE) +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(57, 57, 57, 83, 63, 78, 103, 72, 95, 89), label = "p.signif", method = "t.test", size = 2.5) +
  labs(x = "Population", y = "SROH (Mb)", title ="A") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 35, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black")) 
#plot(plotAll1)

plotA1 <- ggplot(data = data, mapping=aes(x=pop_order, y=A/(1e6), colour=pop_order, fill=pop_order)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) + 
  geom_point(show.legend = FALSE) +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(57, 57, 57, 83, 63, 78, 103, 72, 95, 89), label = "p.signif", method = "t.test", size = 2.5) +
  labs(x = "Population", y = "SROH (Mb)", title ="B") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 35, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black")) 
#plot(plotAll1)

plotB1 <- ggplot(data = data, mapping=aes(x=pop_order, y=B/(1e6), colour=pop_order, fill=pop_order)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) + 
  geom_point(show.legend = FALSE) +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(57, 57, 57, 83, 63, 78, 103, 72, 95, 89), label = "p.signif", method = "t.test", size = 2.5) +
  labs(x = "Population", y = "SROH (Mb)", title ="C") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 35, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black")) 
#plot(plotAll1)

plotC1 <- ggplot(data = data, mapping=aes(x=pop_order, y=C/(1e6), colour=pop_order, fill=pop_order)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE) + 
  geom_point(show.legend = FALSE) +
  #stat_compare_means(comparisons = my_comparisons, label.y = c(57, 57, 57, 83, 63, 78, 103, 72, 95, 89), label = "p.signif", method = "t.test", size = 2.5) +
  labs(x = "Population", y = "SROH (Mb)", title ="D") +
  scale_colour_manual(values = leg) +
  scale_fill_manual(values = alpha(leg, .3)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 35, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black")) 
#plot(plotAll1)

#Plotting and saving
ROHlength <- grid.arrange(plotAll1, plotA1, plotB1, plotC1, ncol=2)
```

![](roh_distribution_5pop_paper_final_files/figure-markdown_github/plotSROH-1.png)

``` r
#ggsave("Total_length_ROH.png", plot = ROHlength, width = 25, height = 29, units = "cm")

#cairo_ps(file = "Total_length_ROH.eps", onefile = FALSE, width=9.8, height=11.4)
#plot(ROHlength)
#dev.off()
```

##### 2. NROH vs SROH (Figure 1)

``` r
#Setting legend colors so it matches Liu. et al (2017) paper
leg <- c("mulatta" = "chartreuse3","lasiotis" = "orchid3", "brevicaudus" = "darkorange", "littoralis" = "red2", "tcheliensis" = "blue3")

ggplot(data) +
  geom_point(mappin=aes(x=TOTAL/(1e6),y=nTOTAL,colour=pop), size = 2.5) + 
  labs(x = "SROH (in Mb)", y = "NROH", title ="", fill="Population") +
  scale_colour_manual(values = leg, name="Subspecies") +
  theme_light() +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), axis.text = element_text(size = 14, color = "black"), axis.ticks = element_line(color = "black"), plot.title = element_blank(), legend.position = c(.82, .23), legend.title = element_text(size=16), legend.text = element_text(size=14,face = 'italic'), legend.key.height = unit(0.4, 'cm'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  coord_cartesian(ylim=c(0,max(data$nTOTAL)), xlim=c(0,max(data$TOTAL/(1e6))))
```

![](roh_distribution_5pop_paper_final_files/figure-markdown_github/srohnrohplot-1.png)

``` r
#ggsave("NROH_vs_SROH.png", width = 15, height = 11, units = "cm")
#dev.off()
#cairo_ps(file = "NROH_vs_SROH.eps", onefile = FALSE, width=5.9, height=4.3)
```

##### 3. Individual ROH genome coverage (Figure 2)

Computing ROH coverage.

``` r
options( "digits"=3)
pop <- c("mulatta","lasiotis", "brevicaudus", "littoralis", "tcheliensis")
genome_size=2900000000
data2 <- cbind(rownames(data),data)
colnames(data2) <- c("ind",colnames(data))
data2 <- data2[,c(1,2,4,6,8,10)]
data2[["nonroh"]] <- genome_size-data$TOTAL
data2[["GiA"]] <- data$A/genome_size*100
data2[["GiB"]] <- data$B/genome_size*100
data2[["GiC"]] <- data$C/genome_size*100
data2[["GiR"]] <- data$TOTAL/genome_size*100
rownames(data2) <- NULL

for (p in pop) {
  df <- data2[data2$pop==p,]
  df <- arrange(df, desc(nonroh))
  assign(paste("data2_", p,sep=""),df)
}

data3 <- rbind(data2_mulatta, data2_lasiotis, data2_littoralis, data2_brevicaudus, data2_tcheliensis)
data3$num <- c(1:79)

longer_data <- data3 %>% #formatting data frame  
        pivot_longer(c(A,B,C,nonroh,GiA,GiB,GiC,GiR), names_to = "class", values_to = "length")
longer_data$ind <- factor(longer_data$ind, levels = data3$ind)

kable(data3[1:10,c(1,2,8,9,10,11)], col.names = c("Individual","Population", 'Gi,a' ,'Gi,b','Gi,c','Gi,r'))
```

| Individual | Population |  Gi,a |  Gi,b |  Gi,c | Gi,r |
|:-----------|:-----------|------:|------:|------:|-----:|
| C_rhe_78   | mulatta    | 0.822 | 0.780 | 0.275 | 1.88 |
| C_rhe_75   | mulatta    | 0.933 | 0.967 | 0.188 | 2.09 |
| C_rhe_72   | mulatta    | 0.885 | 0.871 | 0.355 | 2.11 |
| C_rhe_79   | mulatta    | 0.862 | 0.848 | 0.521 | 2.23 |
| C_rhe_74   | mulatta    | 0.869 | 0.964 | 0.481 | 2.31 |
| C_rhe_70   | mulatta    | 0.857 | 0.905 | 0.666 | 2.43 |
| C_rhe_76   | mulatta    | 0.922 | 0.881 | 0.787 | 2.59 |
| C_rhe_77   | mulatta    | 0.976 | 0.914 | 1.006 | 2.90 |
| C_rhe_73   | mulatta    | 0.981 | 1.024 | 0.997 | 3.00 |
| C_rhe_71   | mulatta    | 1.032 | 1.111 | 1.263 | 3.41 |

``` r
legclass <- c("A" = "steelblue4","B" = "steelblue3", "C" = "steelblue2", "nonroh" = "gray80")
legpop = ifelse(data3$pop %in% "mulatta", "chartreuse3", ifelse(data3$pop %in% "lasiotis", "orchid3",  ifelse(data3$pop %in% "littoralis", "red2",  ifelse(data3$pop %in% "brevicaudus", "darkorange",  "blue3"))))
x_tick <- c(5.5,26,55.5,72,77)

# ROH COVERAGE
rohcoverage.plot <- ggplot(longer_data) +
  geom_col(aes(x=num,y=(length*100/genome_size),fill=class, group = pop)) + 
  geom_segment(aes(y = -1, yend=-1, x = 0.55, xend = 79.5), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 0.55, xend = 0.55), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 10.5, xend = 10.5), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 41.5, xend = 41.5), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 69.5, xend = 69.5), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 74.5, xend = 74.5), size=0.2) +
  geom_segment(aes(y = -0.5, yend=-1, x = 79.5, xend = 79.5), size=0.2) +
  labs(x = "Individuals", y = "Genome coverage (%)", title ="Genome ROH coverage") +
  scale_fill_manual(values = legclass, labels = c("Short ROH","Medium ROH","Long ROH", "Non ROH"), name = "Size class ROH") +
  scale_x_continuous(breaks = x_tick, labels = c("mulatta","lasiotis","littoralis","brevicaudus","tcheliensis")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 30, hjust=0.95, vjust=0.2, size=10, colour = "black", margin = unit(c(-0.9, 0, 0, 0), "cm"), face = 'italic'), axis.text.y = element_text(size = 10, color = "black"), axis.ticks.y = element_line(color = "black"),  axis.title.y = element_text(size = 12), plot.title = element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(),  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.text = element_text(size=8), legend.position = c(.169, .875), legend.title = element_blank(), legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.5, 'cm'), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.line.x = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0.2, 0.2, 1, 0.2), "cm")) +
  coord_cartesian(ylim=c(0,30))

#ggsave("ROH_gen_coverage.png", width = 15, height = 9, units = "cm")
#cairo_ps(file = "ROH_gen_coverage.eps", onefile = FALSE, width=5.5, height=3.54)
#dev.off()
rohcoverage.plot
```

![](roh_distribution_5pop_paper_final_files/figure-markdown_github/save%20plot%20for%20paper-1.png)
