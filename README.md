## shapealign [![Build Status](https://travis-ci.org/thodorisgkolias/shapealign.svg?branch=master)](https://travis-ci.org/thodorisgkolias/shapealign)   [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

### Documentation

'shapealign' is a R package that performs protein alignment using statistical shape analysis techiniques.

### Installation

'shapealign' depends on the Bioconductor packages of 'Biostrings' and 'msa' which can be installed using
```{.r}
source('https://bioconductor.org/biocLite.R')
biocLite('Biostrings')
biocLite('msa')
```

Next, 'shapealign' can be installed from the github repository using 
```{.r}
install.packages('devtools')
devtools::install_github('thodorisgkolias/shapealign')
```

and it can be loaded in R as usually with
```{.r}
library(shapealign)
```


### Code example
Load data 
```{.r}
data1 <- LoadPDB('2gb1')
data2 <- LoadPDB('1ubq')
```
then putting it in a list using

```{.r}
data <- list(data1,data2)
```
Finally you can produce an alignment of these two protein molecules using

```{.r}
ProteinAlign(data, SP = 7, n.cores = 8)
```

### Authors

### License
GPL (>=2)

