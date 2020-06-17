HDplot - Paralog identification using population-level NGS data
================
17 June, 2020

<!-- README.md is generated from README.Rmd. Please edit that file -->

HDplot allows the identification of paralogous SNPs by taking advantage
of the fact that paralogs and non-paralogs differ in expected
heterozygosity and read ratios within heterozygous individuals.

## Input Data

HDplot needs takes a vcfR oject as input.

## Running HDplot

The easiest way to generate the correct input is to load a vcf file
using the read.vcfR() function from the vcfR package.

``` r
library(vcfR)
```

    ## 
    ##    *****       ***   vcfR   ***       *****
    ##    This is vcfR 1.10.0 
    ##      browseVignettes('vcfR') # Documentation
    ##      citation('vcfR') # Citation
    ##    *****       *****      *****       *****

``` r
source("HDplot.R")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
vcfInput<-read.vcfR("exampleInput_100Samples.vcf")
```

    ## Scanning file to determine attributes.
    ## File attributes:
    ##   meta lines: 9
    ##   header_line: 10
    ##   variant count: 19432
    ##   column count: 109
    ## Meta line 9 read in.
    ## All meta lines processed.
    ## gt matrix initialized.
    ## Character matrix gt created.
    ##   Character matrix gt rows: 19432
    ##   Character matrix gt cols: 109
    ##   skip: 0
    ##   nrows: 19432
    ##   row_num: 0
    ## Processed variant 1000Processed variant 2000Processed variant 3000Processed variant 4000Processed variant 5000Processed variant 6000Processed variant 7000Processed variant 8000Processed variant 9000Processed variant 10000Processed variant 11000Processed variant 12000Processed variant 13000Processed variant 14000Processed variant 15000Processed variant 16000Processed variant 17000Processed variant 18000Processed variant 19000Processed variant: 19432
    ## All variants processed

``` r
HDplotResults<-HDplot(vcfInput)
```

The resulting dataframe includes all necessary data to plot results and
assign paralog status

``` r
head(HDplotResults)
```

    ##   CHROM POS ID depth_a depth_b     ratio num_hets num_samples num_called H_all
    ## 1    un  65 11      61      63 0.4919355        4         100         97  0.04
    ## 2    un  17 14      45      32 0.5844156        3         100         98  0.03
    ## 3    un  69 15    1779    2824 0.3864871       93         100         93  0.93
    ## 4    un  73 17     103     104 0.4975845       12         100         95  0.12
    ## 5    un  93 20     319     269 0.5425170       16         100         96  0.16
    ## 6    un  13 21      73      87 0.4562500        7         100         97  0.07
    ##            H       std           D
    ## 1 0.04123711  5.567764  -0.1796053
    ## 2 0.03061224  4.387482   1.4814875
    ## 3 1.00000000 33.922706 -15.4026626
    ## 4 0.12631579  7.193747  -0.0695048
    ## 5 0.16666667 12.124356   2.0619652
    ## 6 0.07216495  6.324555  -1.1067972

The following columns are output in the dataframe: `CHROM, POS, ID`:
These are the chromosome, position, and ID columns from the VCF file.
`depth_a, depth_b`: These are the total read depths for the A and B
alleles in heterozygous individuals. `ratio`: This is the ratio fo the A
reads to B reads. `num_hets, num_samples, num_called`: These are the
number of heterozygous individuals, the number of total individuals, and
the number of individuals with called genotypes. `H_all, H`: H\_all is
the heterozygosity calculated as the number of heterozygous individuals
out of the total individuals. H is the heterozygosity calculated as the
number of heterozygous individuals out of the number of individuals with
called genotypes. With low missing data these values should be similar,
at high missing data H\_all will be reduced relative to H. `std`: This
is the standard deviation of the ratio. `D`: This is the deviation of
the ratio from 50:50 expectations, calculated as a z-score.

## Plotting HDplot results

The important fields for plotting are H and D (ratio can also be an
informative alternative to D). H is the heterozygosity for each SNP
while D is the deviation from even read-ratios in heterozygous
individuals. The ratio column is read-ratio in heterzygous individuals
and is useful for confirming that thresholds set for D look appropriate.

``` r
#plot H and D
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=D))
```

    ## Warning: Removed 127 rows containing missing values (geom_point).

![](tools/Plotting%20Results-1.png)<!-- -->

``` r
#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))
```

    ## Warning: Removed 127 rows containing missing values (geom_point).

![](tools/Plotting%20Results-2.png)<!-- -->
