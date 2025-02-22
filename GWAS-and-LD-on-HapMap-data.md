GWAS and LD on HapMap data
================
Mohanad Hussein
2025-02-16

## OBJECTIVE

This project aims to illustrate an example of how GWAS is manually
performed in R using data from HapMap project. Although this is a part
of the Computational Bio-medicine graduate course at the University of
Göttingen, to me, this projects is a long-waited as I have conducted
web-based GWAS analysis previously, which as expected, did not add much
to the theoretical understanding of Genome Wide Association Studies and
Linkage Disequilibrium. Therefore, I have put some effort in documenting
the analysis by explaining the context behind running GWAS and LD
analysis.

## LIBRARIES AND INPUT DATA

Packages used in this analysis are **snpMatrix** and **snpStats**. The
first is integrates in the second and implements an S4 class to
efficiently handle large SNP data while the latter provides statistical
methods for analyzing SNP data. The tool recommends installing
**hexbin** that enhances plotting options for GWAs convenience and
provides hexagonla bin plots. The data used in this analysis was
retrieved from the HapMap project of which a 1000 samples, varying
between cases (1) and controls (0), where selected for this study. Note
worthy point is that the HapMap project targeted SNPs that have a Minor
Allele Frequency (MAF) \> 0.05 (5%). Each sample is genotyped for 28501
SNPs. The main input is formatted in the form of a matrix where rows are
samples, columns are SNPs and values are numeric codes referring to one
of the three alleles (A,B or AB). The input also requires two additional
elements in addition to the genotype matrix. An element containing
condition information (i.e case/control) and sup-population information
about samples. The other element contains the chromosomal position of
the SNP and the nucleotide base for each of its alleles. Note that this
analysis is dedicated for chromosome 10 only.

``` r
library(snpStats)
library(hexbin)

# loading data, alternatively can be done using load("snps_10.dat")
data(for.exercise)

# inspect the SNP matrix and save the SNP matrix itself
summary(snps.10)
```

    ## $rows
    ##    Call.rate      Certain.calls Heterozygosity  
    ##  Min.   :0.9879   Min.   :1     Min.   :0.0000  
    ##  1st Qu.:0.9896   1st Qu.:1     1st Qu.:0.2993  
    ##  Median :0.9900   Median :1     Median :0.3078  
    ##  Mean   :0.9900   Mean   :1     Mean   :0.3074  
    ##  3rd Qu.:0.9904   3rd Qu.:1     3rd Qu.:0.3159  
    ##  Max.   :0.9919   Max.   :1     Max.   :0.3386  
    ## 
    ## $cols
    ##      Calls        Call.rate     Certain.calls      RAF              MAF        
    ##  Min.   : 975   Min.   :0.975   Min.   :1     Min.   :0.0000   Min.   :0.0000  
    ##  1st Qu.: 988   1st Qu.:0.988   1st Qu.:1     1st Qu.:0.2302   1st Qu.:0.1258  
    ##  Median : 990   Median :0.990   Median :1     Median :0.5030   Median :0.2315  
    ##  Mean   : 990   Mean   :0.990   Mean   :1     Mean   :0.5001   Mean   :0.2424  
    ##  3rd Qu.: 992   3rd Qu.:0.992   3rd Qu.:1     3rd Qu.:0.7671   3rd Qu.:0.3576  
    ##  Max.   :1000   Max.   :1.000   Max.   :1     Max.   :1.0000   Max.   :0.5000  
    ##                                                                                
    ##       P.AA              P.AB             P.BB             z.HWE         
    ##  Min.   :0.00000   Min.   :0.0000   Min.   :0.00000   Min.   :-21.9725  
    ##  1st Qu.:0.06559   1st Qu.:0.2080   1st Qu.:0.06465   1st Qu.: -2.8499  
    ##  Median :0.26876   Median :0.3198   Median :0.27492   Median : -1.1910  
    ##  Mean   :0.34617   Mean   :0.3074   Mean   :0.34647   Mean   : -1.8610  
    ##  3rd Qu.:0.60588   3rd Qu.:0.4219   3rd Qu.:0.60362   3rd Qu.: -0.1014  
    ##  Max.   :1.00000   Max.   :0.5504   Max.   :1.00000   Max.   :  3.7085  
    ##                                                       NA's   :4

``` r
show(snps.10)
```

    ## A SnpMatrix with  1000 rows and  28501 columns
    ## Row names:  jpt.869 ... ceu.464 
    ## Col names:  rs7909677 ... rs12218790

``` r
genotype_matrix = snps.10@.Data
```

## EXPLORATORY DATA ANALYSIS

To have a general overview on the data, we view the head of each
dataframe and summarize the data. The matrix is of a special class
(snpMatrix), so we select only two entries with two SNPs to look at.
After that, we look into the other two input data.

The row and col summary functions compute some summary statistics for
the input matrix. *row.summary* calculates **call rates** and
**heterozygosity** for each sample while *col.names* calculates the
ratio of each of the three alleles per sample in addition to the Z
statistics of Herdy-Weinberg Equilibrium (Z.HWE). The Z.HWE value is a
measue of how far the observed genotype frequency is from the
distribution under HWE principle. A Z.HWE score closer to zero indicates
that observed genotype frequencies are close to what would have been
expected under HWE. Values away from zero (in either direction) marks
significant deviation from HWE due to evolutionary events.

``` r
# inspect the SNP matrix
summary(snps.10)
```

    ## $rows
    ##    Call.rate      Certain.calls Heterozygosity  
    ##  Min.   :0.9879   Min.   :1     Min.   :0.0000  
    ##  1st Qu.:0.9896   1st Qu.:1     1st Qu.:0.2993  
    ##  Median :0.9900   Median :1     Median :0.3078  
    ##  Mean   :0.9900   Mean   :1     Mean   :0.3074  
    ##  3rd Qu.:0.9904   3rd Qu.:1     3rd Qu.:0.3159  
    ##  Max.   :0.9919   Max.   :1     Max.   :0.3386  
    ## 
    ## $cols
    ##      Calls        Call.rate     Certain.calls      RAF              MAF        
    ##  Min.   : 975   Min.   :0.975   Min.   :1     Min.   :0.0000   Min.   :0.0000  
    ##  1st Qu.: 988   1st Qu.:0.988   1st Qu.:1     1st Qu.:0.2302   1st Qu.:0.1258  
    ##  Median : 990   Median :0.990   Median :1     Median :0.5030   Median :0.2315  
    ##  Mean   : 990   Mean   :0.990   Mean   :1     Mean   :0.5001   Mean   :0.2424  
    ##  3rd Qu.: 992   3rd Qu.:0.992   3rd Qu.:1     3rd Qu.:0.7671   3rd Qu.:0.3576  
    ##  Max.   :1000   Max.   :1.000   Max.   :1     Max.   :1.0000   Max.   :0.5000  
    ##                                                                                
    ##       P.AA              P.AB             P.BB             z.HWE         
    ##  Min.   :0.00000   Min.   :0.0000   Min.   :0.00000   Min.   :-21.9725  
    ##  1st Qu.:0.06559   1st Qu.:0.2080   1st Qu.:0.06465   1st Qu.: -2.8499  
    ##  Median :0.26876   Median :0.3198   Median :0.27492   Median : -1.1910  
    ##  Mean   :0.34617   Mean   :0.3074   Mean   :0.34647   Mean   : -1.8610  
    ##  3rd Qu.:0.60588   3rd Qu.:0.4219   3rd Qu.:0.60362   3rd Qu.: -0.1014  
    ##  Max.   :1.00000   Max.   :0.5504   Max.   :1.00000   Max.   :  3.7085  
    ##                                                       NA's   :4

``` r
# inspect SNP support object, info about each snp
head(snp.support)
```

    ##            chromosome position A1 A2
    ## rs7909677          10   101955  A  G
    ## rs7093061          10   112109  C  T
    ## rs12773042         10   117636  C  G
    ## rs7475011          10   133076  C  G
    ## rs11253563         10   148971  A  G
    ## rs4881551          10   149076  A  G

``` r
summary(snp.support)
```

    ##    chromosome    position         A1        A2       
    ##  Min.   :10   Min.   :   101955   A:14019   C: 2349  
    ##  1st Qu.:10   1st Qu.: 28981867   C:12166   G:12254  
    ##  Median :10   Median : 67409719   G: 2316   T:13898  
    ##  Mean   :10   Mean   : 66874497                      
    ##  3rd Qu.:10   3rd Qu.:101966491                      
    ##  Max.   :10   Max.   :135323432

``` r
# inspect subject/sample support object, info about samples, case/ctrl and sub-pop
head(subject.support)
```

    ##         cc stratum
    ## jpt.869  0 JPT+CHB
    ## jpt.862  0 JPT+CHB
    ## jpt.948  0 JPT+CHB
    ## ceu.564  0     CEU
    ## ceu.904  0     CEU
    ## ceu.665  0     CEU

``` r
summary(subject.support)
```

    ##        cc         stratum   
    ##  Min.   :0.0   CEU    :494  
    ##  1st Qu.:0.0   JPT+CHB:506  
    ##  Median :0.5                
    ##  Mean   :0.5                
    ##  3rd Qu.:1.0                
    ##  Max.   :1.0

To understand how the data is distributed, we plot the distribution of
MAF “summarizing” columns in Genotype matrix. MAF is calculated by
dividing the Minor Allele Count by the total alleles at that position
(2N). The second plot show the distribution of Z scores for each SNP
under HWE principal. We see that Z scores are mostly close to zero but
there is a clear deviation from the HWE principal as the distribution is
not symmetric.

``` r
sample_sum = row.summary(snps.10)
snp_sum = col.summary(snps.10)

hist(snp_sum$MAF)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
hist(snp_sum$z.HWE)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

Next, we plot a correlogram to see the correlation between variables of
the summary stat, namely call rates, certain calls and heterozygosity.
An outlier sample appears to have a near-zero value in the
Heteroztygosity column. Since Loss Of Heterozygosity is a sign of
inbreeding and can hide the signal of actually causal variants.

``` r
plot(sample_sum)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

After the correlogram, we plot two other hexagonal binning plots to
check difference of call rates in contrls vs samples, and also check
their MAF differences.

``` r
# compare call rates across condition ctrl/case
cases = subject.support$cc == 1
controls = subject.support$cc == 0

# recompute the genotype column for cases and 
cases_sum = col.summary(snps.10[cases,])
controls_sum = col.summary(snps.10[controls,])

# visualize call rates in hexagonal bin plot
hb = hexbin(controls_sum$Call.rate, cases_sum$Call.rate, xbin = 50) # create the hexbin object
sp = plot(hb) # plot the hexbin object
hexVP.abline(sp$plot.vp, 0, 1, col="black") # add slope
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# visualize allele frequencies in hexagonal bin plot
sp = plot(hexbin(controls_sum$MAF, cases_sum$MAF, xbin = 50))
hexVP.abline(sp$plot.vp, 0,1, col = "white")
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## GWAS ANALYSIS

Some Qc steps are recommended in literature, some of them are
[here](https://pmc.ncbi.nlm.nih.gov/articles/PMC6001694/table/mpr1608-tbl-0001/)

The logic of GWAS is to screen thousands of SNPs to find correlation
with a trait -cases in this example- using Cochran-Armitage test. We
have, for each SNP, chi-squared tests on 1 (additive model of minor
allele effect) and 2 (no assumption about allele) degrees of freedom
(df), together with N, the number of subjects for whom data were
available. The 2 df test is the conventional Pearsonian test for the 3 ×
2 contingency table. The large number of NA values for the latter test
reflects the fact that, for these SNPs, the minor allele frequency was
such that one homozygous genotype did not occur in the data.

Some filtering needs to be done on the GWAS results; we consider taking
only SNPs with MAF value of higher than 1% in the population.
Additionally, we also would like to restrict our analysis to samples
that are under HWE principle and thus filter out those ridiculously out
of HWE.

We filter the GWAS results for the mentioned criteria and store the
positions of the filtered SNPs to use in creating Manhattan plot.

``` r
# run GWAS
tests = single.snp.tests(cc, data = subject.support, snp.data = snps.10)
use2 = snp_sum$MAF > 0.01 & snp_sum$z.HWE^2 < 200
tests = tests[use2]
positions = snp.support[use2, "position"]
```

Another way to look into the significance of SNPs is the QQ plot (not so
relevant here)

So far, this is broad view of SNPs that are in input data that is used
to run GWAS. However, we usually would like to account for population
structure. This can be done by defining a factor of strata in the
subject support object. This requires rerunning GWAS steps and

## LINKAGE DISEQUILIBRIUM

Two genomic positions are said to be in Linkage Disequilibrium (LD) if
their segregation occurred dependent of each other. In other words, they
segregated together. There are different metrics to quantify LD, of
which we choose D.prime and R.squared. The first measures presence of LD
between SNPs (i.e how much do two SNPs appear together in samples) which
infers recombination. The latter measures how strong the correlation
between SNPs is. This measure is illustrated in what is called LD plot;
showing 1 as black and 0 as white.

``` r
data(ld.example)

# run LD with two calculation methods
snps.10_LD = ld(snps.10, stats = c("D.prime", "R.squared"), depth = 100)

# plot LD plots for both metrics
image(snps.10_LD$D.prime[1:500,1:500], lwd = 0)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
image(snps.10_LD$R.squared[1:500, 1:500], lwd = 0)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

## COMBINATION OF GWAS WITH LD ANALYSIS

Each of the two analysis on its own provides comprehension for an
underlying genomic context. However, despite the statistical power of
GWAS, it doesn’t pinpoint the causal SNP/s in a disease. Here comes the
rule of LD; when given a set of SNPs that are significantly correlated
to a trait, we can highlight which of these SNPs fall in LD region and
also quantify it’s strength. A region with high LD is more likely to
contain haplotypes that are passed on to disease carriers. Yet not
enough to definitively determine a causal variant, it provides us with a
relatively narrow genomic range of SNPs to direct for fine mapping. Such
a narrow region can be visaulized as follows.

``` r
# puplishing ready plot of narrow LD region of 200 SNPs selected upon checking higher GWAS hits' positions
spectrum <- rainbow(10, start=0, end=1/6)[10:1]
image(snps.10_LD$D.prime[75:274,75:274], lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE)
```

![](GWAS-and-LD-on-HapMap-data_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

## DATA AVAILABILITY

This analysis is part of the Computational Bio-medicine course at the
University of Göttingen and relies heavily on the documentation provided
by the used R packages. Raw code and data files are available on my
Github: <https://github.com/mohanadhussein/GWAS_LD_analysis>
