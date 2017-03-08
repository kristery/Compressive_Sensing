%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 10000;
interval = 40;
fmax = 800;

t0 = 1:siglen/interval:siglen;
f0 = 1:500/interval:fmax;
sigma = 1:100/interval:100;
        
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
   
%signal = sin(seq).*(cos(log(seq+7)))+seq;
[signal, Fs] = audioread('voice.wav');
signal = signal';

dictionary = real(dictionary);

disp(size(dictionary));
tic
[YFIT,R,COEFF] = wmpalg('OMP', signal, dictionary');
toc

audiowrite('original.wav', YFIT, Fs);
err0 = sum(abs(YFIT-signal'));

tic
[H, W, diff] = NMF_reduce(dictionary, 100);
toc
disp(size(H));
disp(size(W));
%%% 
tic
[YFIT,R,COEFF] = wmpalg('OMP', signal, W');
toc
audiowrite('new.wav', YFIT, Fs);
err1 = sum(abs(YFIT-signal'));
disp(size(COEFF));
disp(err0);
disp(err1);

disp((err1-err0)/err0);

%x = inv(W'*W)*W'*COEFF;
%recovery = (dictionary'-diff)*x;

%disp(COEFF);
%disp(sum(abs(recovery-signal')));
%disp(sum(abs(YFIT-signal)));

