%%% it is used to test matchin pursuit %%%
% parameters: t0, f0, sigma
siglen = 1000;
interval = 40;
fmax = 600;

t0 = 1:siglen/interval:siglen;
f0 = 1:500/interval:fmax;
sigma = 0.1:5/interval:5;
        
dictionary = zeros(length(t0)*length(f0)*length(sigma), siglen);
seq = 1:siglen;
idx0 = 1;
for idx1 = 1:length(t0)
    for idx2 = 1:length(f0)
        for idx3 = 1:length(sigma)
            dictionary(idx0,:) = (2^0.25)/(sigma(idx3)^0.5)*exp(j*2*pi*f0(idx2)*seq-pi*((seq-t0(idx1)).^2)/(sigma(idx3)^2));
            idx0 = idx0+1;
        end
    end
end
   
signal = sin(seq).*(cos(log(seq)));

dictionary = real(dictionary);

tic
[YFIT,R,COEFF] = wmpalg('BMP', signal, dictionary');
toc

%disp(COEFF);
disp(size(dictionary));
disp(size(signal));
disp(sum(abs(YFIT)));