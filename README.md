# snpReady
A tool to assist breeders to prepare genotypic datasets for genomic analysis in order to run genomic analysis and estimates some population genetics parameters. Thus, it produce outputs that can be use in many packages or softwares related to genomic analysis.

## Installation

snpReady is available on CRAN
```R
install.packages("snpReady")
```
In github it is available the experimental version and its installation needs to be done via [devtools](https://github.com/hadley/devtools#updating-to-the-latest-version-of-devtools). Hence, it is necessary first install devtools and later install snpReady
```R
install.packages("devtools")
library(devtools)
install_github("italo-granato/snpReady")
```

## Usage
Below, we present some basic usage for the three functions available in snpReady

### raw.data

Function to clean and recode raw dataset from genotyping

```R
data(maize.line)
M <- raw.data(as.matrix(maize.line), frame="long", base=TRUE, sweep.sample= 0.8, 
call.rate=0.95, maf=0.05, input=TRUE, outfile="-101")

```
### G.matrix

Function to create genomic relationship matrix (GRM)

```R
data(maize.hyb)
x <- G.matrix(maize.hyb, method = "VanRaden", format = "wide")
A <- x$Ga
D <- x$Gd
```
### popgen

Function to estimate some parameters of genetic of population using markers

```R
data(maize.hyb)
x <- popgen(maize.hyb) 
```
 

## Acknowledgments

I would like to thank people from [Allogamous Plant Breeding Laboratory Team](http://www.genetica.esalq.usp.br/alogamas/index2.html) for helping 
in this project. 

## Authors

[Allogamous Plant Breeding Laboratory Team](http://www.genetica.esalq.usp.br/alogamas/index2.html)

### Contributting
To anyone who wants to contribute, please contact [Italo Granato](mailto:italo.granato@gmail.com) for more details.
