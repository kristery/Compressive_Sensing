a = rand(300,500);

opt = statset('UseParallel', true);

tic
[W,H] = nnmf(a,200);
disp(size(W));
disp(size(H));
toc
tic
[W,H] =  nnmf(a, 200, 'options', opt);
toc

