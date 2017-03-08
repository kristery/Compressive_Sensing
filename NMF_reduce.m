%%% Use NMF to lower the dimension so that MP could be feasible %%%
function [H, W, diff] = NMF_reduce(A, c)
    diff = min(min(A));
    pA = A - diff;
	k = c;
	n = size(A, 1);
	m = size(A, 2);
	Winit = rand(n, k);
	Hinit = rand(k, m);
    max_iter = 50;
    %%% Find H, W such that min(|A-HW|)
    addpath('./NMF_CD/sparse_square');
	%[H, W] = nnmf(pA,c);
	[H W] = sparse_CD(pA, c, max_iter, Winit, Hinit, 1);
end
