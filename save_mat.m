%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 1000;
interval = 97;

tmax = 10000;
dt = tmax/interval;
t = 1:dt:tmax;

fmax = 4000;
df = fmax/97;
f = 1:df:fmax;

smax = 300;
ds = 300/7;
sigma = 1:ds:smax;

dictionary = zeros(length(t)*2*length(f)*length(sigma), siglen);
n = 1:siglen;
idx0 = 1;
for idx1 = 1:length(t)
	for idx2 = 1:length(f)
		for idx3 = 1:length(sigma)
			dictionary(2*idx0-1,:) = (2^0.25)/(sigma(idx3)^0.5)*cos(2*pi*f(idx2)*n)./exp(pi*((n-t(idx1)).^2)/(sigma(idx3)^2));
			dictionary(2*idx0,:) = (2^0.25)/(sigma(idx3)^0.5)*sin(2*pi*f(idx2)*n)./exp(pi*((n-t(idx1)).^2)/(sigma(idx3)^2));
			idx0 = idx0+1;
		end
	end
end

disp(size(dictionary));
dictionary = sparse(dictionary);
disp('Begin to find unique set of dictionary.')
dictionary = unique(dictionary,'rows');
disp(size(dictionary));
save('unique_dict/dict1000_97.mat', 'dictionary', '-v7.3');

disp('Save dictionary successful!')

[H, W, diff] = NMF_reduce(dictionary, 300);
H = sparse(H);
W = sparse(W);
save('unique_dict/H1000_97_300.mat','H', '-v7.3');
save('unique_dict/W1000_97_300.mat','W', '-v7.3');
