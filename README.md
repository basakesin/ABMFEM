# Annealing-Based Model-Free Expectation Maximization Algorithm (ABMFEM)

 this code provides an example run of the anneling-based model-free expectation maximization algorithm provided by the
 abmfem package. 
 
 If you use this code, please cite:

 Köktürk, Başak Esin, and Bilge Karaçalı. "Annealing-based model-free expectation maximisation for
 multi-colour flow cytometry data clustering." International Journal of Data Mining and Bioinformatics 14,
 no. 1 (2016): 86-99.

 This package have reqired the qsl toolbox to estimate the posterior
 probabilities. qsl package is avaliable on: http://web.iyte.edu.tr/~bilgekaracali/Projects/QSL/qsl.tar.gz


## Usage
----
```
library(mvtnorm)
n<-1000
x<- rmvnorm(n, c(0,0))
y<- rmvnorm(n, c(0,4))
z <- rbind(x,y)
label1 <-rep(0, times = n)
label2 <-rep(1, times = n)
label<- c(label1,label2)
D <- as.matrix(dist(z, method = "euclidean", p = 2, diag = TRUE, upper = TRUE))
N <- 2*n
k<- 1
P0<-qsl(z,label,D,N,k)
```
## Dependencies
-----
In order to use qsl function, you need to add kstar.R function into your workspace.

## Author
----
Basak Esin KOKTURK GUZEL

## Credits
----
QSL toolbox for MATLAB

http://likya.iyte.edu.tr/~bilgekaracali/Projects/QSL/


If you use the algorithm, please cite:

B. Karacali, "Quasi-supervised learning for biomedical data analysis," Pattern Recognition, Vol. 43, No. 10, p. 3674-3682, 2010
