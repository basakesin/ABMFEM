# Annealing-Based Model-Free Expectation Maximization Algorithm (ABMFEM) MATLAB Package

 This code provides an example run of the anneling-based model-free expectation maximization algorithm provided by the
 abmfem package. 
 
Annealing-Based Model-Free Expectation Maximization Algorithm can be estimate the number of clusters and assigned samples to this clusters without any pre-knowledge. 

## Usage
----
```octave
X0 = randn(1000,2);
X1 = [randn(1000,2)+repmat([0 4],1000,1)];
X2 = [randn(1000,2)+repmat([4 4],1000,1)]

figure(1)
clf
plot(X0(:,1),X0(:,2),'.');
hold on
plot(X1(:,1),X1(:,2),'r.');
plot(X2(:,1),X2(:,2),'g.')
xlabel('first variate');
ylabel('second variate');
legend('x in C_0','x in C_1','x in C_2');
title('Vectors in C_0, C_1 and C_2');
axis equal
grid on

X=[X0; X1; X2];
Y=[zeros(1000,1); ones(1000,1); 2*ones(1000,1)];

[data, labels, clusters] = divideData(X);

figure(2)
gscatter(data(:,1),data(:,2),labels)
xlabel('first variate');
ylabel('second variate');
legend('x in C_0','x in C_1','x in C_2');
title('Estimated Clusters using ABMFEM Algorithm');
axis equal
grid on

```
## Dependencies
-----
 This package have reqired the qsl toolbox to estimate the posterior
 
 probabilities. qsl package is avaliable on: http://web.iyte.edu.tr/~bilgekaracali/Projects/QSL/qsl.tar.gz

## Author
----
Basak Esin KOKTURK GUZEL

## Credits
----

 If you use this code, please cite:

 Köktürk, Başak Esin, and Bilge Karaçalı. "Annealing-based model-free expectation maximisation for
 multi-colour flow cytometry data clustering." International Journal of Data Mining and Bioinformatics 14,
 no. 1 (2016): 86-99.
 
## License
----
 The abmfem package is free software; you can redistribute it and/or  modify it under the terms of the GNU General Public License as 
 published by the Free Software Foundation; either version 3 of the License, or any later version. 
 
The software package is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
