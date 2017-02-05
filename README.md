## shapealign

### Documentation

'shapealign' is a R package that performs protein alignment using statistical shape analysis techiniques.

### Installation
The package can be install from the github repository using 
```{.r}
install.packages('devtools')
devtools::install_github('thodorisgkolias/shapealign')
```

The packages of 'Biostrings' and 'msa' can be installed using
```{.r}
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("msa")
```
then it can be loaded in R as usually with
```{.r}
library(shapealign)
```



### Code example
Load data 
```{.r}
data1 = LoadPDB('2gb1')
data2 = LoadPDB('1ubq')
```
then putting it in a list using

```{.r}
data = list(data1,data2)
```
Finally you can produce an alignment of these two protein molecules using

```{.r}
ProteinAlign(data, SP = 7, n.cores = 8)
```

### Authors

### License
GPL-2

