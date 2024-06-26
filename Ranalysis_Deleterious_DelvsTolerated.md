### **Deleterious vs tolerated variation**

#### 1. Importing the genotypes for each individual

``` r
#Import genotypes
delsites_genotypes <- read.table("deleterious_genotypes_5pop2.txt", header=TRUE, sep=c("\t"),stringsAsFactors = TRUE)
tolsites_genotypes <- read.table("tolerated_genotypes_5pop2.txt", header=TRUE, sep=c("\t"),stringsAsFactors = TRUE)

#Import population file
population_file <- read.table("79Chinese.id2pop.2.txt", header=FALSE, sep="\t",stringsAsFactors = TRUE)
colnames(population_file) <- c("ind","pop")
population_file$pop <- factor(population_file$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
population_file <- population_file %>% arrange(pop)
population_file <- population_file %>%
  mutate(ind = factor(ind, ind))

genotypes_df <- rbind.data.frame(setNames(cbind.data.frame(rep("deleterious",nrow(delsites_genotypes)),delsites_genotypes),c("pred",colnames(delsites_genotypes))), setNames(cbind.data.frame(rep("tolerated", nrow(tolsites_genotypes)), tolsites_genotypes),c("pred",colnames(tolsites_genotypes))))

kable(genotypes_df[1:10, 1:10])
```

| pred        | CHROM |       POS | C_rhe_6 | C_rhe_7 | C_rhe_8 | C_rhe_9 | C_rhe_10 | C_rhe_11 | C_rhe_12 |
|:---------|:-----|--------:|:------|:------|:------|:------|:-------|:-------|:-------|
| deleterious | chr1  |   1218571 | 1/1     | 1/1     | 1/1     | ./.     | 1/1      | 1/1      | ./.      |
| deleterious | chr1  |  85200868 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      |
| deleterious | chr1  | 114893039 | 0/1     | 0/1     | 0/1     | 0/1     | 0/1      | 0/1      | 0/0      |
| deleterious | chr2  |  11017841 | 0/0     | 0/0     | ./.     | 0/0     | 0/0      | 0/0      | 0/0      |
| deleterious | chr2  |  15900832 | 0/0     | 0/0     | 0/0     | 0/0     | ./.      | ./.      | ./.      |
| deleterious | chr2  |  27572545 | ./.     | 0/0     | 0/0     | 0/0     | 0/0      | 0/1      | 0/0      |
| deleterious | chr2  |  27572871 | 0/0     | 0/0     | 0/0     | 0/0     | 0/0      | 0/0      | 0/0      |
| deleterious | chr2  |  27746394 | 0/0     | 0/0     | 0/0     | 0/0     | 0/1      | ./.      | ./.      |
| deleterious | chr2  |  50956358 | 0/0     | 0/0     | ./.     | 0/0     | 0/0      | 0/0      | 0/0      |
| deleterious | chr2  |  51007414 | 0/0     | 0/1     | 0/0     | 0/0     | 0/0      | 0/0      | 0/1      |

#### 2. Counting the number of tolerated and deleterious genotypes (0/0, 0/1, 1/1) for each individual

``` r
#counting the genotypes
genotypes_count_del <- setNames(cbind.data.frame(rep("deleterious",79), as.data.frame(colnames(delsites_genotypes)[3:81])),c("pred","ind"))
genotypes_count_del$del.homoref <- colSums(delsites_genotypes[,3:ncol(delsites_genotypes)] == "0/0", na.rm = TRUE)
genotypes_count_del$del.hetero <- colSums(delsites_genotypes[,3:ncol(delsites_genotypes)] == "0/1", na.rm = TRUE)
genotypes_count_del$del.homoalt <- colSums(delsites_genotypes[,3:ncol(delsites_genotypes)] == "1/1", na.rm = TRUE)
genotypes_count_del$del.miss <- colSums(delsites_genotypes[,3:ncol(delsites_genotypes)] == "./.", na.rm = TRUE)

genotypes_count_tol <- setNames(cbind.data.frame(rep("tolerated",79), as.data.frame(colnames(tolsites_genotypes)[3:81])),c("pred","ind"))
genotypes_count_tol$tol.homoref <- colSums(tolsites_genotypes[,3:ncol(tolsites_genotypes)] == "0/0", na.rm = TRUE)
genotypes_count_tol$tol.hetero <- colSums(tolsites_genotypes[,3:ncol(tolsites_genotypes)] == "0/1", na.rm = TRUE)
genotypes_count_tol$tol.homoalt <- colSums(tolsites_genotypes[,3:ncol(tolsites_genotypes)] == "1/1", na.rm = TRUE)
genotypes_count_tol$tol.miss <- colSums(tolsites_genotypes[,3:ncol(tolsites_genotypes)] == "./.", na.rm = TRUE)

#adding the population name
genotypes_count_del$pop <- NA
for (i in 1:nrow(genotypes_count_del)) {
  ind=as.character(genotypes_count_del$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  genotypes_count_del[genotypes_count_del$ind==ind, "pop"] <- pop
}
genotypes_count_tol$pop <- NA
for (i in 1:nrow(genotypes_count_tol)) {
  ind=as.character(genotypes_count_tol$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  genotypes_count_tol[genotypes_count_tol$ind==ind, "pop"] <- pop
}

genotypes_count_del <- as.data.frame(unclass(genotypes_count_del), stringsAsFactors = TRUE)
genotypes_count_tol <- as.data.frame(unclass(genotypes_count_tol), stringsAsFactors = TRUE)
#with(genotypes_count_del, genotypes_count_del[order(pop),])
#with(genotypes_count_tol, genotypes_count_tol[order(pop),])

genotypes_count <- cbind.data.frame(genotypes_count_del[,2],genotypes_count_del[,7], genotypes_count_del[,c(3,4,5,6)], genotypes_count_tol[,c(3,4,5,6)])
colnames(genotypes_count) <- c("ind","pop","del.homoref","del.hetero","del.homoalt", "del.miss", "tol.homoref","tol.hetero","tol.homoalt", "tol.miss")
genotypes_count <- genotypes_count %>%
    mutate(total = select(., del.homoref:tol.miss) %>% rowSums(na.rm = TRUE))

genotypes_count$pop <- factor(genotypes_count$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
genotypes_count <- genotypes_count %>% arrange(pop)
genotypes_count <- genotypes_count %>%
  mutate(pop, ind = factor(ind, ind))

#print
#write.csv(genotypes_count, file = "genotypes_counts.csv")
kable(genotypes_count[1:10,], col.names = c("Individual","Population", "del.0/0","del.0/1","del.1/1","del./.","tol.0/0","tol.0/1","tol.1/1","tol./.","total"))
```

| Individual | Population | del.0/0 | del.0/1 | del.1/1 | del./. | tol.0/0 | tol.0/1 | tol.1/1 | tol./. | total |
|:--------|:--------|------:|------:|------:|-----:|------:|------:|------:|-----:|-----:|
| C_rhe_70   | mulatta    |    1931 |     200 |      33 |    221 |   10937 |    1953 |     583 |   1440 | 17298 |
| C_rhe_71   | mulatta    |    2073 |     209 |      42 |     61 |   11952 |    1894 |     653 |    414 | 17298 |
| C_rhe_74   | mulatta    |    2063 |     245 |      37 |     40 |   12036 |    2036 |     628 |    213 | 17298 |
| C_rhe_75   | mulatta    |    2034 |     231 |      34 |     86 |   11702 |    2053 |     659 |    499 | 17298 |
| C_rhe_79   | mulatta    |    2063 |     251 |      26 |     45 |   11988 |    1980 |     639 |    306 | 17298 |
| C_rhe_78   | mulatta    |    2039 |     256 |      39 |     51 |   11929 |    2011 |     642 |    331 | 17298 |
| C_rhe_76   | mulatta    |    1999 |     228 |      37 |    121 |   11559 |    2030 |     664 |    660 | 17298 |
| C_rhe_77   | mulatta    |    2062 |     234 |      30 |     59 |   11808 |    2032 |     653 |    420 | 17298 |
| C_rhe_73   | mulatta    |    2081 |     227 |      39 |     38 |   11904 |    2014 |     670 |    325 | 17298 |
| C_rhe_72   | mulatta    |    2057 |     239 |      34 |     55 |   11922 |    2009 |     632 |    350 | 17298 |

*Figure 3*

``` r
leg <- c("mulatta" = "chartreuse3","lasiotis" = "orchid3", "brevicaudus" = "darkorange", "littoralis" = "red2", "tcheliensis" = "blue3")

#DELETERIOUS HOMOZYGOTES
delhetero <- ggplot(genotypes_count) +
  geom_violin(mapping = aes(x=pop, y=del.homoalt, fill = pop, colour = pop), trim=TRUE, scale="area", show.legend = FALSE) + 
  geom_point(aes(x=pop, y=del.homoalt, fill = pop, colour = pop), show.legend = FALSE) +
  scale_y_continuous(breaks = seq(20, 70, 5), limits = c(20, 70)) +
  labs(x = "Individuals", y = "Deleterious homozygotes", title ="A") +
  scale_fill_manual(values = alpha(leg, 0.3)) +
  scale_colour_manual(values = leg) +
  guides(fill=guide_legend(title="Subspecies")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

#DELETERIOUS HETEROZYGOTES
delhomo <- ggplot(genotypes_count) +
  geom_violin(mapping = aes(x=pop, y=del.hetero, fill = pop, colour = pop), trim=TRUE, scale="area", show.legend = FALSE) + 
  geom_point(aes(x=pop, y=del.hetero, fill = pop, colour = pop), show.legend = FALSE) +
  labs(x = "Individuals", y = "Deleterious heterozygotes", title ="B") +
  scale_y_continuous(breaks = seq(100, 500, 25), limits = c(150, 300)) +
  scale_fill_manual(values = alpha(leg, 0.3)) +
  scale_colour_manual(values = leg) +
  guides(fill=guide_legend(title="Subspecies")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

#DELETERIOUS ALLELES
delallele <- ggplot(genotypes_count) +
  geom_violin(mapping = aes(x=pop, y=del.hetero+2*del.homoalt, fill = pop, colour = pop), trim=TRUE, scale="area", show.legend = FALSE) + 
  geom_point(aes(x=pop, y=del.hetero+2*del.homoalt, fill = pop, colour = pop), show.legend = FALSE) +
  labs(x = "Individuals", y = "Deleterious alles", title ="C") +
  scale_y_continuous(breaks = seq(100, 500, 25), limits = c(225, 400)) +
  scale_fill_manual(values = alpha(leg, 0.3)) +
  scale_colour_manual(values = leg) +
  guides(fill=guide_legend(title="Subspecies")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12, angle = 40, hjust = 1, face = "italic"), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

TOTDEL <- grid.arrange(delhetero, delhomo, delallele, ncol=3)
```

<img src="Deleterious_distribution-_in_ROH_paper_files/figure-markdown_github/heterozhistogramprop-1.png" style="display: block; margin: auto;" />

``` r
#ggsave("Total_deleterious_variation.png", plot = TOTDEL, width = 30, height = 11, units = "cm")
#cairo_ps(file = "Total_deleterious_variation.eps", onefile = FALSE, width=13, height=4.33)
#plot(TOTDEL)
#dev.off()
```

Wilcox-test (two-sided or alternative = less), deleterious homozygotes

``` r
wtest_homoalt <- c()
for (pop1 in c("mulatta", "lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
  pop1
  for (pop2 in c("mulatta","lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
    x <- genotypes_count[genotypes_count$pop == pop1, ]$del.homoalt
    y <- genotypes_count[genotypes_count$pop == pop2, ]$del.homoalt
    t2s <- wilcox.test(x = x, y = y, paired = FALSE)
    tl <- wilcox.test(x = x, y = y, paired = FALSE, alternative = "less")
    wtest_homoalt <- c(wtest_homoalt, pop1, pop2, formatC(t2s$p.value, digits = 5), formatC(tl$p.value, digits = 5))
  }
}
wtest_homoalt_mat <- matrix(data = wtest_homoalt, ncol = 4, byrow = TRUE)
colnames(wtest_homoalt_mat) <- c("pop1", "pop2", "p, two-sided", "p, less")
wtest_homoalt_mat
```

    ##       pop1          pop2          p, two-sided p, less     
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51527"   
    ##  [2,] "mulatta"     "lasiotis"    "0.00084655" "0.00042328"
    ##  [3,] "mulatta"     "littoralis"  "0.015326"   "0.0076632" 
    ##  [4,] "mulatta"     "brevicaudus" "0.016635"   "0.0083175" 
    ##  [5,] "mulatta"     "tcheliensis" "0.0082837"  "0.0041418" 
    ##  [6,] "lasiotis"    "mulatta"     "0.00084655" "0.99962"   
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50283"   
    ##  [8,] "lasiotis"    "littoralis"  "0.64795"    "0.68148"   
    ##  [9,] "lasiotis"    "brevicaudus" "0.040993"   "0.020497"  
    ## [10,] "lasiotis"    "tcheliensis" "0.14793"    "0.073966"  
    ## [11,] "littoralis"  "mulatta"     "0.015326"   "0.99301"   
    ## [12,] "littoralis"  "lasiotis"    "0.64795"    "0.32397"   
    ## [13,] "littoralis"  "littoralis"  "     1"     "0.50328"   
    ## [14,] "littoralis"  "brevicaudus" "0.032529"   "0.016265"  
    ## [15,] "littoralis"  "tcheliensis" "0.11841"    "0.059207"  
    ## [16,] "brevicaudus" "mulatta"     "0.016635"   "0.99409"   
    ## [17,] "brevicaudus" "lasiotis"    "0.040993"   "0.98167"   
    ## [18,] "brevicaudus" "littoralis"  "0.032529"   "0.98567"   
    ## [19,] "brevicaudus" "brevicaudus" "     1"     "0.54224"   
    ## [20,] "brevicaudus" "tcheliensis" "0.52836"    "0.79974"   
    ## [21,] "tcheliensis" "mulatta"     "0.0082837"  "0.99714"   
    ## [22,] "tcheliensis" "lasiotis"    "0.14793"    "0.93226"   
    ## [23,] "tcheliensis" "littoralis"  "0.11841"    "0.9465"    
    ## [24,] "tcheliensis" "brevicaudus" "0.52836"    "0.26418"   
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54224"

Wilcox-test (two-sided or alternative = less), deleterious heterozygotes

``` r
wtest_hetero <- c()
for (pop1 in c("mulatta", "lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
  pop1
  for (pop2 in c("mulatta", "lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
    x <- genotypes_count[genotypes_count$pop == pop1, ]$del.hetero
    y <- genotypes_count[genotypes_count$pop == pop2, ]$del.hetero
    t2s <- wilcox.test(x = x, y = y, paired = FALSE)
    tl <- wilcox.test(x = x, y = y, paired = FALSE, alternative = "less")
    wtest_hetero <- c(wtest_hetero, pop1, pop2, formatC(t2s$p.value, digits = 5), formatC(tl$p.value, digits = 5))
  }
}
wtest_hetero_mat <- matrix(data = wtest_hetero, ncol = 4, byrow = TRUE)
colnames(wtest_hetero_mat) <- c("pop1", "pop2", "p, two-sided", "p, less")
wtest_hetero_mat
```

    ##       pop1          pop2          p, two-sided p, less     
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51513"   
    ##  [2,] "mulatta"     "lasiotis"    "0.77277"    "0.62519"   
    ##  [3,] "mulatta"     "littoralis"  "0.61881"    "0.70219"   
    ##  [4,] "mulatta"     "brevicaudus" "0.90244"    "0.59684"   
    ##  [5,] "mulatta"     "tcheliensis" "0.002664"   "0.99933"   
    ##  [6,] "lasiotis"    "mulatta"     "0.77277"    "0.38639"   
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50282"   
    ##  [8,] "lasiotis"    "littoralis"  "0.8315"     "0.59017"   
    ##  [9,] "lasiotis"    "brevicaudus" "0.78337"    "0.6258"    
    ## [10,] "lasiotis"    "tcheliensis" "0.00058886" "0.99975"   
    ## [11,] "littoralis"  "mulatta"     "0.61881"    "0.3094"    
    ## [12,] "littoralis"  "lasiotis"    "0.8315"     "0.41575"   
    ## [13,] "littoralis"  "littoralis"  "     1"     "0.50327"   
    ## [14,] "littoralis"  "brevicaudus" "0.90005"    "0.56979"   
    ## [15,] "littoralis"  "tcheliensis" "0.0010015"  "0.99958"   
    ## [16,] "brevicaudus" "mulatta"     "0.90244"    "0.45122"   
    ## [17,] "brevicaudus" "lasiotis"    "0.78337"    "0.39168"   
    ## [18,] "brevicaudus" "littoralis"  "0.90005"    "0.45003"   
    ## [19,] "brevicaudus" "brevicaudus" "     1"     "0.54224"   
    ## [20,] "brevicaudus" "tcheliensis" "0.015873"   "0.99603"   
    ## [21,] "tcheliensis" "mulatta"     "0.002664"   "0.001332"  
    ## [22,] "tcheliensis" "lasiotis"    "0.00058886" "0.00029443"
    ## [23,] "tcheliensis" "littoralis"  "0.0010015"  "0.00050075"
    ## [24,] "tcheliensis" "brevicaudus" "0.015873"   "0.0079365" 
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54224"

Wilcox-test (two-sided or alternative = less), deleterious alleles

``` r
wtest_allel <- c()
for (pop1 in c("mulatta", "lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
  for (pop2 in c("mulatta", "lasiotis", "littoralis", "brevicaudus", "tcheliensis")) {
    x <- genotypes_count[genotypes_count$pop == pop1, ]$del.hetero + 2*genotypes_count[genotypes_count$pop == pop1, ]$del.homoalt
    y <- genotypes_count[genotypes_count$pop == pop2, ]$del.hetero + 2*genotypes_count[genotypes_count$pop == pop2, ]$del.homoalt
    t2s <- wilcox.test(x = x, y = y, paired = FALSE)
    tl <- wilcox.test(x = x, y = y, paired = FALSE, alternative = "less")
    wtest_allel <- c(wtest_allel, pop1, pop2, formatC(t2s$p.value, digits = 5), formatC(tl$p.value, digits = 5))
  }
}
wtest_allel_mat <- matrix(data = wtest_allel, ncol = 4, byrow = TRUE)
colnames(wtest_allel_mat) <- c("pop1", "pop2", "p, two-sided", "p, less")
wtest_allel_mat
```

    ##       pop1          pop2          p, two-sided p, less    
    ##  [1,] "mulatta"     "mulatta"     "     1"     "0.51513"  
    ##  [2,] "mulatta"     "lasiotis"    "0.0056513"  "0.0028256"
    ##  [3,] "mulatta"     "littoralis"  "0.056542"   "0.028271" 
    ##  [4,] "mulatta"     "brevicaudus" "0.007992"   "0.003996" 
    ##  [5,] "mulatta"     "tcheliensis" "0.51349"    "0.78022"  
    ##  [6,] "lasiotis"    "mulatta"     "0.0056513"  "0.99743"  
    ##  [7,] "lasiotis"    "lasiotis"    "     1"     "0.50282"  
    ##  [8,] "lasiotis"    "littoralis"  "0.32329"    "0.84205"  
    ##  [9,] "lasiotis"    "brevicaudus" "0.036929"   "0.018465" 
    ## [10,] "lasiotis"    "tcheliensis" "0.073719"   "0.96669"  
    ## [11,] "littoralis"  "mulatta"     "0.056542"   "0.97381"  
    ## [12,] "littoralis"  "lasiotis"    "0.32329"    "0.16165"  
    ## [13,] "littoralis"  "littoralis"  "     1"     "0.50327"  
    ## [14,] "littoralis"  "brevicaudus" "0.039484"   "0.019742" 
    ## [15,] "littoralis"  "tcheliensis" "0.13849"    "0.9372"   
    ## [16,] "brevicaudus" "mulatta"     "0.007992"   "0.99767"  
    ## [17,] "brevicaudus" "lasiotis"    "0.036929"   "0.98351"  
    ## [18,] "brevicaudus" "littoralis"  "0.039484"   "0.98254"  
    ## [19,] "brevicaudus" "brevicaudus" "     1"     "0.54224"  
    ## [20,] "brevicaudus" "tcheliensis" "0.095238"   "0.97222"  
    ## [21,] "tcheliensis" "mulatta"     "0.51349"    "0.25674"  
    ## [22,] "tcheliensis" "lasiotis"    "0.073719"   "0.03686"  
    ## [23,] "tcheliensis" "littoralis"  "0.13849"    "0.069246" 
    ## [24,] "tcheliensis" "brevicaudus" "0.095238"   "0.047619" 
    ## [25,] "tcheliensis" "tcheliensis" "     1"     "0.54224"

#### 3. Deleterious vs tolerated variants in ROH

``` r
del_distrib <- read.table("deleterious_distribution_5pop.txt", header = TRUE, stringsAsFactors = TRUE)
tol_distrib <- read.table("tolerated_distribution_5pop.txt", header = TRUE, stringsAsFactors = TRUE)

options(digits=3)

pop <- c("mulatta","lasiotis", "brevicaudus", "littoralis", "tcheliensis")

#Import population file
population_file <- read.table("79Chinese.id2pop.2.txt", header=FALSE, sep="\t",stringsAsFactors = TRUE)
colnames(population_file) <- c("ind","pop")
population_file$pop <- factor(population_file$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
population_file <- population_file %>% arrange(pop)
population_file <- population_file %>%
  mutate(ind = factor(ind, ind))

#Add population to individuals
del_distrib$pop <- NA
for (i in 1:nrow(del_distrib)) {
  ind=as.character(del_distrib$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  del_distrib[del_distrib$ind==ind, "pop"] <- pop
}
tol_distrib$pop <- NA
for (i in 1:nrow(tol_distrib)) {
  ind=as.character(tol_distrib$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  tol_distrib[tol_distrib$ind==ind, "pop"] <- pop
}

#Reorder individuals
del_distrib$pop <- factor(del_distrib$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
del_distrib <- del_distrib %>% arrange(pop)
del_distrib <- del_distrib %>%
  mutate(pop, ind = factor(ind, ind))
tol_distrib$pop <- factor(tol_distrib$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
tol_distrib <- tol_distrib %>% arrange(pop)
tol_distrib <- tol_distrib %>%
  mutate(pop, ind = factor(ind, ind))

#Creating summary for del and tol
distrib_summary <- rbind(setNames(cbind(rep("deleterious", 79), del_distrib),c("pred",colnames(del_distrib))), setNames(cbind(rep("tolerated", 79), tol_distrib),c("pred",colnames(tol_distrib))))

distrib_summary[["deleterious"]] <- as.numeric(distrib_summary$pred=='deleterious')
distrib_summary[["frac.inA"]] <- distrib_summary$nA/distrib_summary$ntotal
distrib_summary[["frac.inB"]] <- distrib_summary$nB/distrib_summary$ntotal
distrib_summary[["frac.inC"]] <- distrib_summary$nC/distrib_summary$ntotal
distrib_summary[["frac.inany"]] <- distrib_summary$nallroh/distrib_summary$ntotal

options(digits=3)
#write.csv(distrib_summary, file="DelTol_distribinROH.csv")
kable(distrib_summary[1:10,])
```

| pred        | ind      |  nA |  nB |  nC | nallroh | ntotal | pop     | deleterious | frac.inA | frac.inB | frac.inC | frac.inany |
|:-------|:------|--:|--:|--:|-----:|----:|:-----|-------:|------:|------:|------:|-------:|
| deleterious | C_rhe_70 |   1 |   0 |   0 |       1 |     33 | mulatta |           1 |    0.030 |    0.000 |    0.000 |      0.030 |
| deleterious | C_rhe_71 |   2 |   2 |   0 |       4 |     42 | mulatta |           1 |    0.048 |    0.048 |    0.000 |      0.095 |
| deleterious | C_rhe_72 |   1 |   2 |   1 |       4 |     34 | mulatta |           1 |    0.029 |    0.059 |    0.029 |      0.118 |
| deleterious | C_rhe_73 |   2 |   1 |   1 |       4 |     39 | mulatta |           1 |    0.051 |    0.026 |    0.026 |      0.103 |
| deleterious | C_rhe_74 |   0 |   2 |   0 |       2 |     37 | mulatta |           1 |    0.000 |    0.054 |    0.000 |      0.054 |
| deleterious | C_rhe_75 |   2 |   0 |   0 |       2 |     34 | mulatta |           1 |    0.059 |    0.000 |    0.000 |      0.059 |
| deleterious | C_rhe_76 |   1 |   3 |   1 |       5 |     37 | mulatta |           1 |    0.027 |    0.081 |    0.027 |      0.135 |
| deleterious | C_rhe_77 |   0 |   0 |   0 |       0 |     30 | mulatta |           1 |    0.000 |    0.000 |    0.000 |      0.000 |
| deleterious | C_rhe_78 |   0 |   0 |   0 |       0 |     39 | mulatta |           1 |    0.000 |    0.000 |    0.000 |      0.000 |
| deleterious | C_rhe_79 |   1 |   1 |   0 |       2 |     26 | mulatta |           1 |    0.038 |    0.038 |    0.000 |      0.077 |

*Figure S3*

``` r
distrib_summary_trans <- distrib_summary[, c(1,2,3,4,5,6,7,8)] %>%
  pivot_longer(cols = c(nA, nB, nC, nallroh), names_to = "class", values_to = "n")
distrib_summary_trans$class <- factor(distrib_summary_trans$class, levels = c("nallroh","nA", "nB","nC"))

#DELETERIOUS HOMO IN ROH
delinroh <- ggplot(distrib_summary_trans[distrib_summary_trans$pred=="deleterious",]) +
  geom_boxplot(mapping = aes(x=class, y=n, fill = pop, colour = pop)) + 
  #geom_point(aes(x=class, y=n, color = pop, group = interaction(class, pop)), position = position_dodge(width = 0.75), show.legend = FALSE) + 
  labs(x = "Population", y = "Deleterious homozygotes in ROH", title ="A") +
  scale_y_continuous(breaks = seq(0, 50, 5), limits = c(0, 50)) +
  scale_fill_manual(values = alpha(leg, 0.5)) +
  scale_colour_manual(values = leg) +
  guides(fill=guide_legend(title="Subspecies", nrow = 1), colour=guide_legend(title="Subspecies", nrow = 1)) +
  scale_x_discrete(labels = c("Any size ROH", "Short ROH", "Medium ROH", "Long ROH")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.28, .90), legend.title = element_text(size=12), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 

#delinroh
#ggsave("Deleterious_variation_in_ROH.png", plot = delinroh, width = 28, height = 10, units = "cm")
#cairo_ps(file = "Deleterious_variation_in_ROH.eps", onefile = FALSE, width=11.02, height=3.94)
#plot(delinroh)
#dev.off()

#TOLERATED HOMO IN ROH
tolinroh <- ggplot(distrib_summary_trans[distrib_summary_trans$pred=="tolerated",]) +
  geom_boxplot(mapping = aes(x=class, y=n, fill = pop, colour = pop), show.legend = FALSE) + 
  #geom_point(aes(x=class, y=n, color = pop, group = interaction(class, pop)), position = position_dodge(width = 0.75), show.legend = FALSE) + 
  labs(x = "Population", y = "Tolerated homozygotes in ROH", title ="B") +
  scale_y_continuous(breaks = seq(0, 500, 50), limits = c(0, 400)) +
  scale_fill_manual(values = alpha(leg, 0.5)) +
  scale_colour_manual(values = leg) +
  guides(fill=guide_legend(title="Subspecies")) +
  scale_x_discrete(labels = c("Any size ROH", "Short ROH", "Medium ROH", "Long ROH")) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 12), axis.title.x=element_blank(), axis.text.y = element_text(size=10, color = "black"), axis.title.y=element_text(size=14, margin = unit(c(0, 0.5, 0,0), "cm")), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.12, .85), legend.title = element_text(size=14), legend.text = element_text(size=11), legend.key.height = unit(0.35, 'cm'), plot.title = element_text(size= 15), plot.margin = unit(c(0.2, 1, 0.2, 0.75), "cm"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), panel.border = element_blank()) 

#tolinroh
#ggsave("Tolerated_variation_in_ROH.png", plot = tolinroh, width = 28, height = 10, units = "cm")
#cairo_ps(file = "Tolerated_variation_in_ROH.eps", onefile = FALSE, width=11.02, height=3.94)
#plot(tolinroh)
#dev.off()
```

#### 4. Importing ROH coverage per individual

``` r
roh <- read.table("5pop_roh_summary.txt", header=TRUE, sep="\ ", stringsAsFactors = TRUE)
roh[["ind"]] <- rownames(roh)
rownames(roh) <- NULL
genome_size=2900000000

#Add population names
roh$pop <- NA
for (i in 1:nrow(roh)) {
  ind=as.character(roh$ind[i])
  pop=as.character(population_file[population_file$ind==ind, "pop"])
  roh[roh$ind==ind, "pop"] <- pop
}

#Reorder individuals
roh$pop <- factor(roh$pop, levels = c("mulatta","lasiotis","brevicaudus","littoralis","tcheliensis"))
roh <- roh %>% arrange(pop)
roh <- roh %>%
  mutate(pop, ind = factor(ind, ind))

distrib_summary[["cov.any"]] <- NA 
distrib_summary[["cov.A"]] <- NA
distrib_summary[["cov.B"]] <- NA
distrib_summary[["cov.C"]] <- NA

for (i in 1:nrow(distrib_summary)) {
  ind=distrib_summary[i,]$ind
  distrib_summary[i,"cov.any"] <- roh[roh$ind==ind,][1,"TOTAL"]/genome_size
  distrib_summary[i,"cov.A"] <- roh[roh$ind==ind,][1,"A"]/genome_size
  distrib_summary[i,"cov.B"] <- roh[roh$ind==ind,][1,"B"]/genome_size
  distrib_summary[i,"cov.C"] <- roh[roh$ind==ind,][1,"C"]/genome_size
}

kable(distrib_summary[1:10,], col.names = c("Prediction","Population","Individual","Class A","Class B"," Class C", "Any class" , "Total", "Deleterious", "Fraction in A","Fraction in B"," Fraction in C", "Fraction in any class","Any ROH coverage","A coverage", "B coverage","C coverage"))
```

| Prediction  | Population | Individual | Class A | Class B | Class C | Any class | Total   | Deleterious | Fraction in A | Fraction in B | Fraction in C | Fraction in any class | Any ROH coverage | A coverage | B coverage | C coverage |
|:----|:---|---:|---:|---:|---:|---:|:---|----:|----:|----:|----:|------:|-----:|---:|---:|---:|
| deleterious | C_rhe_70   |          1 |       0 |       0 |       1 |        33 | mulatta |           1 |         0.030 |         0.000 |         0.000 |                 0.030 |            0.024 |      0.009 |      0.009 |      0.007 |
| deleterious | C_rhe_71   |          2 |       2 |       0 |       4 |        42 | mulatta |           1 |         0.048 |         0.048 |         0.000 |                 0.095 |            0.034 |      0.010 |      0.011 |      0.013 |
| deleterious | C_rhe_72   |          1 |       2 |       1 |       4 |        34 | mulatta |           1 |         0.029 |         0.059 |         0.029 |                 0.118 |            0.021 |      0.009 |      0.009 |      0.004 |
| deleterious | C_rhe_73   |          2 |       1 |       1 |       4 |        39 | mulatta |           1 |         0.051 |         0.026 |         0.026 |                 0.103 |            0.030 |      0.010 |      0.010 |      0.010 |
| deleterious | C_rhe_74   |          0 |       2 |       0 |       2 |        37 | mulatta |           1 |         0.000 |         0.054 |         0.000 |                 0.054 |            0.023 |      0.009 |      0.010 |      0.005 |
| deleterious | C_rhe_75   |          2 |       0 |       0 |       2 |        34 | mulatta |           1 |         0.059 |         0.000 |         0.000 |                 0.059 |            0.021 |      0.009 |      0.010 |      0.002 |
| deleterious | C_rhe_76   |          1 |       3 |       1 |       5 |        37 | mulatta |           1 |         0.027 |         0.081 |         0.027 |                 0.135 |            0.026 |      0.009 |      0.009 |      0.008 |
| deleterious | C_rhe_77   |          0 |       0 |       0 |       0 |        30 | mulatta |           1 |         0.000 |         0.000 |         0.000 |                 0.000 |            0.029 |      0.010 |      0.009 |      0.010 |
| deleterious | C_rhe_78   |          0 |       0 |       0 |       0 |        39 | mulatta |           1 |         0.000 |         0.000 |         0.000 |                 0.000 |            0.019 |      0.008 |      0.008 |      0.003 |
| deleterious | C_rhe_79   |          1 |       1 |       0 |       2 |        26 | mulatta |           1 |         0.038 |         0.038 |         0.000 |                 0.077 |            0.022 |      0.009 |      0.008 |      0.005 |

*Figure 4*

``` r
leg <- c("deleterious" = "red","tolerated" = "blue")

#Any class
plot_list <- list()
plot_list[[1]] <- ggplot_gtable(ggplot_build(ggplot() +
  geom_point(data=distrib_summary, mapping = aes(x=cov.any, y=frac.inany, colour = pred, shape = pop)) +
  guides(shape = "none") +
  scale_colour_manual(values = leg, name="SIFT Prediction", labels = c("Deleterious","Tolerated")) +
  geom_smooth(data=distrib_summary[distrib_summary$pred=="deleterious",], mapping = aes(x=cov.any, y=frac.inany), colour="red", method=lm, se = FALSE, size = 0.5) +
  geom_smooth(data=distrib_summary[distrib_summary$pred=="tolerated",], mapping = aes(x=cov.any, y=frac.inany), colour="blue", method=lm, se = FALSE, size = 0.5) +
  labs(x = "Fraction of genome in any size ROH", y = paste(expression("Fraction of alternate\nhomozygotes in any size ")," ROH", sep=""), colour="SIFT Prediction", shape="Population", title = "A") +
  theme_light() +
  theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size = 10, margin = unit(c(0, 0.4, 0, 0), "cm")), axis.text.x = element_text(size=8, color="black"), axis.ticks = element_line(color="black"), axis.text.y = element_text(size = 8, color="black"), legend.position = c(.22, .85), legend.title = element_text(size=9), legend.text = element_text(size=8, margin = unit(c(0, 0.2, 0,0), "cm")), legend.key.height = unit(0.3, 'cm'), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))))

plot <- ggplot() +
    geom_point(data=distrib_summary, mapping = aes(x=cov.A, y=frac.inA, colour = pred, shape = pop), show.legend = FALSE) + 
    scale_colour_manual(values = leg, name="SIFT Prediction") +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="deleterious",], mapping = aes(x=cov.A, y=frac.inA), colour="red", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="tolerated",], mapping = aes(x=cov.A, y=frac.inA), colour="blue", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    labs(x = "Fraction of genome in short ROH", y = "Fraction of alternate\nhomozygotes in short ROH", colour="SIFT Prediction", shape="Population", title = "B") +
    theme_light() +
    theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size = 10, margin = unit(c(0, 0.4, 0, 0), "cm")), axis.text.x = element_text(size=8, color="black"), axis.ticks = element_line(color="black"), axis.text.y = element_text(size = 8, color="black"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 
gplot <- ggplot_gtable(ggplot_build(plot))
plot_list[[2]] <- gplot  

plot <- ggplot() +
    geom_point(data=distrib_summary, mapping = aes(x=cov.B, y=frac.inB, colour = pred, shape = pop), show.legend = FALSE) + 
    scale_colour_manual(values = leg, name="SIFT Prediction") +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="deleterious",], mapping = aes(x=cov.B, y=frac.inB), colour="red", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="tolerated",], mapping = aes(x=cov.B, y=frac.inB), colour="blue", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    labs(x = "Fraction of genome in medium ROH", y = "Fraction of alternate\nhomozygotes in medium ROH", colour="SIFT Prediction", shape="Population", title = "C") +
    theme_light() +
    theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size = 10, margin = unit(c(0, 0.4, 0, 0), "cm")), axis.text.x = element_text(size=8, color="black"), axis.ticks = element_line(color="black"), axis.text.y = element_text(size = 8, color="black"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 
gplot <- ggplot_gtable(ggplot_build(plot))
plot_list[[3]] <- gplot  

plot <- ggplot() +
    geom_point(data=distrib_summary, mapping = aes(x=cov.C, y=frac.inC, colour = pred, shape = pop)) + 
    scale_colour_manual(values = leg, name="SIFT Prediction") +
   guides(colour = "none") +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="deleterious",], mapping = aes(x=cov.C, y=frac.inC), colour="red", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    geom_smooth(data=distrib_summary[distrib_summary$pred=="tolerated",], mapping = aes(x=cov.C, y=frac.inC), colour="blue", method=lm, se = FALSE, size = 0.5, show.legend = FALSE) +
    labs(x ="Fraction of genome in long ROH", y = "Fraction of alternate\nhomozygotes in medium ROH", colour="SIFT Prediction", shape="Population", title = "D") +
    theme_light() +
    theme(axis.title.x = element_text(size=10), axis.title.y = element_text(size = 10, margin = unit(c(0, 0.4, 0, 0), "cm")), axis.text.x = element_text(size=8, color="black"), axis.text.y = element_text(size = 8, color="black"), axis.ticks = element_line(color="black"), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), plot.title = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_blank(), legend.position = c(.22, .78), legend.title = element_text(size=9), legend.text = element_text(size=8, margin = unit(c(0, 0.2, 0,0), "cm")), legend.key.height = unit(0.3, 'cm'), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) 
gplot <- ggplot_gtable(ggplot_build(plot))
plot_list[[4]] <- gplot

fig <- do.call("grid.arrange", c(plot_list, ncol=2))
```

<img src="Deleterious_distribution-_in_ROH_paper_files/figure-markdown_github/plotinoutROH-1.png" style="display: block; margin: auto;" />

``` r
#ggsave("Frac_DEL_ROH.png", fig, width = 20, height = 16, units = "cm")

#cairo_ps(file = "Frac_DEL_ROH.eps", onefile = FALSE, width=7.87, height=6.3)
#plot(fig)
#dev.off()
```

**Pearson correlations**

``` r
pearson_cor <- setNames(as.data.frame(matrix(ncol = 4, nrow = 8)), c("ROH class", "Prediction", "Pearson r", "p-value")) 
pearson_cor[1,] <- c("AnyROH", "Deleterious", cor(distrib_summary[distrib_summary$pred=="deleterious",]$cov.any, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inany, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="deleterious",]$cov.any, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inany, method = "pearson")$p.value)
pearson_cor[2,] <- c("AnyROH", "Tolerated", cor(distrib_summary[distrib_summary$pred=="tolerated",]$cov.any, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inany, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="tolerated",]$cov.any, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inany, method = "pearson")$p.value)
pearson_cor[3,] <- c("A", "Deleterious", cor(distrib_summary[distrib_summary$pred=="deleterious",]$cov.A, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inA, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="deleterious",]$cov.A, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inA, method = "pearson")$p.value)
pearson_cor[4,] <- c("A", "Tolerated", cor(distrib_summary[distrib_summary$pred=="tolerated",]$cov.A, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inA, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="tolerated",]$cov.A, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inA, method = "pearson")$p.value)
pearson_cor[5,] <- c("B", "Deleterious", cor(distrib_summary[distrib_summary$pred=="deleterious",]$cov.B, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inB, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="deleterious",]$cov.B, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inB, method = "pearson")$p.value)
pearson_cor[6,] <- c("B", "Tolerated", cor(distrib_summary[distrib_summary$pred=="tolerated",]$cov.B, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inB, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="tolerated",]$cov.B, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inB, method = "pearson")$p.value)
pearson_cor[7,] <- c("C", "Deleterious", cor(distrib_summary[distrib_summary$pred=="deleterious",]$cov.C, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inC, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="deleterious",]$cov.C, distrib_summary[distrib_summary$pred=="deleterious",]$frac.inC, method = "pearson")$p.value)
pearson_cor[8,] <- c("C", "Tolerated", cor(distrib_summary[distrib_summary$pred=="tolerated",]$cov.C, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inC, method = "pearson"), cor.test(distrib_summary[distrib_summary$pred=="tolerated",]$cov.C, distrib_summary[distrib_summary$pred=="tolerated",]$frac.inC, method = "pearson")$p.value)

options(digits = 3)
kable(pearson_cor, digits=3)
```

| ROH class | Prediction  | Pearson r         | p-value              |
|:----------|:------------|:------------------|:---------------------|
| AnyROH    | Deleterious | 0.877837290233746 | 2.4843831241627e-26  |
| AnyROH    | Tolerated   | 0.961939518766159 | 3.89358952607375e-45 |
| A         | Deleterious | 0.235658780925201 | 0.0365515010003264   |
| A         | Tolerated   | 0.540395538281519 | 2.74001737934822e-07 |
| B         | Deleterious | 0.627705827077362 | 5.97356700131289e-10 |
| B         | Tolerated   | 0.868246239747514 | 3.78402629605971e-25 |
| C         | Deleterious | 0.90663434593545  | 1.3841087081563e-30  |
| C         | Tolerated   | 0.957194911425552 | 3.28297778466173e-43 |

**Linear regressions**

*Any size ROH*

``` r
mod.any <- lm(data=distrib_summary, frac.inany ~ cov.any*deleterious)
summary(mod.any)
```

    ## 
    ## Call:
    ## lm(formula = frac.inany ~ cov.any * deleterious, data = distrib_summary)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14537 -0.03150 -0.00403  0.02898  0.22645 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)          0.04382    0.01164    3.76  0.00024 ***
    ## cov.any              1.83677    0.11945   15.38  < 2e-16 ***
    ## deleterious         -0.00836    0.01647   -0.51  0.61225    
    ## cov.any:deleterious  0.70608    0.16892    4.18  4.9e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0571 on 154 degrees of freedom
    ## Multiple R-squared:  0.824,  Adjusted R-squared:  0.82 
    ## F-statistic:  240 on 3 and 154 DF,  p-value: <2e-16

*Class A ROH*

``` r
modA <- lm(data=distrib_summary, frac.inA ~ cov.A*deleterious) 
summary(modA)
```

    ## 
    ## Call:
    ## lm(formula = frac.inA ~ cov.A * deleterious, data = distrib_summary)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.05667 -0.01176 -0.00206  0.00966  0.09787 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        0.01486    0.01071    1.39    0.167  
    ## cov.A              2.01143    0.93070    2.16    0.032 *
    ## deleterious        0.00107    0.01515    0.07    0.944  
    ## cov.A:deleterious  0.68431    1.31621    0.52    0.604  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0278 on 154 degrees of freedom
    ## Multiple R-squared:  0.0986, Adjusted R-squared:  0.0811 
    ## F-statistic: 5.62 on 3 and 154 DF,  p-value: 0.00111

*Class B ROH*

``` r
modB <- lm(data=distrib_summary, frac.inB ~ cov.B*deleterious) 
summary(modB)
```

    ## 
    ## Call:
    ## lm(formula = frac.inB ~ cov.B * deleterious, data = distrib_summary)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.07228 -0.01275 -0.00080  0.00922  0.08804 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.01011    0.00607    1.67    0.098 .  
    ## cov.B              2.36661    0.31584    7.49  4.9e-12 ***
    ## deleterious        0.00665    0.00858    0.78    0.439    
    ## cov.B:deleterious  0.59974    0.44666    1.34    0.181    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0272 on 154 degrees of freedom
    ## Multiple R-squared:  0.508,  Adjusted R-squared:  0.499 
    ## F-statistic:   53 on 3 and 154 DF,  p-value: <2e-16

*Class C ROH*

``` r
modC <- lm(data=distrib_summary, frac.inC ~ cov.C*deleterious)
summary(modC)
```

    ## 
    ## Call:
    ## lm(formula = frac.inC ~ cov.C * deleterious, data = distrib_summary)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.16558 -0.01780 -0.00338  0.01741  0.16167 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.00641    0.00745    0.86     0.39    
    ## cov.C              1.86904    0.10415   17.95  < 2e-16 ***
    ## deleterious       -0.00992    0.01053   -0.94     0.35    
    ## cov.C:deleterious  0.62870    0.14729    4.27  3.4e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0437 on 154 degrees of freedom
    ## Multiple R-squared:  0.855,  Adjusted R-squared:  0.852 
    ## F-statistic:  303 on 3 and 154 DF,  p-value: <2e-16
