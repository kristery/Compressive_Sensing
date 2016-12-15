%%% Use NMF to lower the dimension so that MP could be feasible %%%
function [H, W, diff] = NMF_reduce(A, c)
    diff = min(min(A));
    pA = A - diff;
    
	opt = statset('useParallel', true);    
    %%% Find H, W such that min(|A-HW|)
    [H, W] = nnmf(pA,c,'options', opt);
end
