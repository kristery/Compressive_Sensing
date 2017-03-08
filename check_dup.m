%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 10000;

% A = HW, A' = W'H'
dictionary = load('unique_dict/dict_formal1.mat');
dictionary = dictionary.dictionary;
dictionary = real(dictionary);

%H = load('unique_dict/H_formal1_100.mat');
%H = H.H;

%W = load('unique_dict/W_formal1_100.mat');
%W = W.W;
disp(size(dictionary));
disp(size(unique(dictionary,'rows')));
