N = 1000;
Na = 300;
A = matrix(runif(9*N^2),3*N, 3*N);
Ei = matrix(runif(3*N*Na),3*N, Na);
system.time(E <-  solve(A, Ei))
# user  system elapsed 
# 18.155   0.095  18.395 
