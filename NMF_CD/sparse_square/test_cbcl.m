load ../data/cbcl;

%V = sparse(V);
V = sprand(1000,10000,0.2);

n = size(V,1);
m = size(V,2);
disp(size(V));
k = 100;

Winit = rand(n,k);
Hinit = rand(k,m);
maxiter = 20;
type = 1; %% GCD
tic
[W H] = sparse_CD(V, k, maxiter, Winit, Hinit, type);
toc
sum(sum(abs(V-W*H)))

tic
[W, H] = nnmf(V, k);
toc
sum(sum(abs(V-W*H)))
