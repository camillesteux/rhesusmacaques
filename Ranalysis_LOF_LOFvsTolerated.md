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
kable(STOP_genotypes[1:10, 1:10])
```

| CHROM |      POS | C_rhe_6 | C_rhe_7 | C_rhe_8 | C_rhe_9 | C_rhe_10 | C_rhe_11 | C_rhe_12 | C_rhe_13 |
|:-----|-------:|:------|:------|:------|:------|:-------|:-------|:-------|:-------|
| chr2  | 27745889 | 0/0     | 0/0     | 0/0     | 0/0     | ./.      | 0/0      | 0/0      | ./.      |
| chr2  | 50713433 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 50713636 | 0/0     | 0/0     | 1/1     | 0/1     | 0/1      | 1/1      | 0/0      | 0/1      |
| chr2  | 50713683 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 50713927 | 0/1     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 57956388 | 0/0     | 0/0     | 0/1     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |
| chr2  | 84368555 | 0/0     | ./.     | 0/0     | 0/0     | 1/1      | ./.      | 0/0      | ./.      |
| chr3  |   541039 | 0/0     | 0/0     | ./.     | 0/0     | 0/0      | ./.      | 0/0      | 0/0      |
| chr3  | 33007093 | ./.     | 0/0     | ./.     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |
| chr3  | 40003876 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      | 0/0      |

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

Wilcox-test (two-sided or alternative = less), LOF homozygotes

    ##       pop1          pop2          p, two-sided p, less   
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51587" 
    ##  [2,] "mulatta"     "lasiotis"    "0.44473"    "0.22236" 
    ##  [3,] "mulatta"     "littoralis"  "0.76265"    "0.6314"  
    ##  [4,] "mulatta"     "brevicaudus" "0.80121"    "0.4006"  
    ##  [5,] "mulatta"     "tcheliensis" "0.23651"    "0.90463" 
    ##  [6,] "lasiotis"    "mulatta"     "0.44473"    "0.78682" 
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50288" 
    ##  [8,] "lasiotis"    "littoralis"  "0.26237"    "0.87206" 
    ##  [9,] "lasiotis"    "brevicaudus" "0.83343"    "0.60142" 
    ## [10,] "lasiotis"    "tcheliensis" "0.14283"    "0.93472" 
    ## [11,] "littoralis"  "mulatta"     "0.76265"    "0.38132" 
    ## [12,] "littoralis"  "lasiotis"    "0.26237"    "0.13118" 
    ## [13,] "littoralis"  "littoralis"  "     1"     "0.50331" 
    ## [14,] "littoralis"  "brevicaudus" "0.59342"    "0.29671" 
    ## [15,] "littoralis"  "tcheliensis" "0.59367"    "0.7205"  
    ## [16,] "brevicaudus" "mulatta"     "0.80121"    "0.64717" 
    ## [17,] "brevicaudus" "lasiotis"    "0.83343"    "0.41672" 
    ## [18,] "brevicaudus" "littoralis"  "0.59342"    "0.72064" 
    ## [19,] "brevicaudus" "brevicaudus" "     1"     "0.54333" 
    ## [20,] "brevicaudus" "tcheliensis" "0.39167"    "0.85786" 
    ## [21,] "tcheliensis" "mulatta"     "0.23651"    "0.11825" 
    ## [22,] "tcheliensis" "lasiotis"    "0.14283"    "0.071414"
    ## [23,] "tcheliensis" "littoralis"  "0.59367"    "0.29684" 
    ## [24,] "tcheliensis" "brevicaudus" "0.39167"    "0.19583" 
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54333"

Wilcox-test (two-sided or alternative = less), LOF heterozygotes

    ##       pop1          pop2          p, two-sided p, less    
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51537"  
    ##  [2,] "mulatta"     "lasiotis"    "0.73739"    "0.64274"  
    ##  [3,] "mulatta"     "littoralis"  "0.029051"   "0.98666"  
    ##  [4,] "mulatta"     "brevicaudus" "0.063535"   "0.9761"   
    ##  [5,] "mulatta"     "tcheliensis" "0.026524"   "0.9904"   
    ##  [6,] "lasiotis"    "mulatta"     "0.73739"    "0.3687"   
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50283"  
    ##  [8,] "lasiotis"    "littoralis"  "0.011114"   "0.99468"  
    ##  [9,] "lasiotis"    "brevicaudus" "0.14738"    "0.93252"  
    ## [10,] "lasiotis"    "tcheliensis" "0.028906"   "0.98715"  
    ## [11,] "littoralis"  "mulatta"     "0.029051"   "0.014525" 
    ## [12,] "littoralis"  "lasiotis"    "0.011114"   "0.0055571"
    ## [13,] "littoralis"  "littoralis"  "     1"     "0.50329"  
    ## [14,] "littoralis"  "brevicaudus" "     1"     "0.52016"  
    ## [15,] "littoralis"  "tcheliensis" "0.28784"    "0.86725"  
    ## [16,] "brevicaudus" "mulatta"     "0.063535"   "0.031767" 
    ## [17,] "brevicaudus" "lasiotis"    "0.14738"    "0.07369"  
    ## [18,] "brevicaudus" "littoralis"  "     1"     "   0.5"   
    ## [19,] "brevicaudus" "brevicaudus" "     1"     "0.54224"  
    ## [20,] "brevicaudus" "tcheliensis" "0.46198"    "0.82787"  
    ## [21,] "tcheliensis" "mulatta"     "0.026524"   "0.013262" 
    ## [22,] "tcheliensis" "lasiotis"    "0.028906"   "0.014453" 
    ## [23,] "tcheliensis" "littoralis"  "0.28784"    "0.14392"  
    ## [24,] "tcheliensis" "brevicaudus" "0.46198"    "0.23099"  
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54333"

Wilcox-test (two-sided or alternative = less), LOF alleles

    ##       pop1          pop2          p, two-sided p, less    
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51522"  
    ##  [2,] "mulatta"     "lasiotis"    "0.71393"    "0.35697"  
    ##  [3,] "mulatta"     "brevicaudus" "0.45921"    "0.80607"  
    ##  [4,] "mulatta"     "littoralis"  "0.072701"   "0.96622"  
    ##  [5,] "mulatta"     "tcheliensis" "0.15713"    "0.93795"  
    ##  [6,] "lasiotis"    "mulatta"     "0.71393"    "0.65436"  
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50283"  
    ##  [8,] "lasiotis"    "brevicaudus" "0.26989"    "0.87478"  
    ##  [9,] "lasiotis"    "littoralis"  "0.0054515"  "0.9974"   
    ## [10,] "lasiotis"    "tcheliensis" "0.062609"   "0.9718"   
    ## [11,] "brevicaudus" "mulatta"     "0.45921"    "0.2296"   
    ## [12,] "brevicaudus" "lasiotis"    "0.26989"    "0.13495"  
    ## [13,] "brevicaudus" "brevicaudus" "     1"     "0.54224"  
    ## [14,] "brevicaudus" "littoralis"  "0.66868"    "0.68379"  
    ## [15,] "brevicaudus" "tcheliensis" "0.54762"    "0.78968"  
    ## [16,] "littoralis"  "mulatta"     "0.072701"   "0.03635"  
    ## [17,] "littoralis"  "lasiotis"    "0.0054515"  "0.0027257"
    ## [18,] "littoralis"  "brevicaudus" "0.66868"    "0.33434"  
    ## [19,] "littoralis"  "littoralis"  "     1"     "0.50328"  
    ## [20,] "littoralis"  "tcheliensis" "0.31419"    "0.8547"   
    ## [21,] "tcheliensis" "mulatta"     "0.15713"    "0.078566" 
    ## [22,] "tcheliensis" "lasiotis"    "0.062609"   "0.031305" 
    ## [23,] "tcheliensis" "brevicaudus" "0.54762"    "0.27381"  
    ## [24,] "tcheliensis" "littoralis"  "0.31419"    "0.1571"   
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54224"

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
