%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
rng('shuffle');
siglen = 10000;

% A = HW, A' = W'H'
dictionary = load('unique_dict/dict_formal1.mat');
dictionary = dictionary.dictionary;
dictionary = real(dictionary);

[signal, Fs] = audioread('whale.wav');
signal = signal';

disp(size(dictionary));
Sd = size(dictionary);
tic
[YFIT,R,COEFF, iopt] = wmpalg('BMP', signal, dictionary');
toc

iopt = unique(iopt);

audiowrite('original.wav', YFIT, Fs);
err0 = sum(abs(YFIT-signal'));

tic
idx = randperm(Sd(1),1579);
new_dict = dictionary(idx,:);
disp(size(new_dict));
[YFIT0,R,COEFF,iopt1] = wmpalg('BMP', signal, new_dict');
toc

iopt1 = unique(iopt1);

disp(size(iopt));
disp(size(iopt1));
disp(size(intersect(iopt,iopt1)));

audiowrite('rand.wav', YFIT0, Fs);

err1 = sum(abs(YFIT0-signal'));

disp(err0);
disp(err1);
disp((err1-err0)/err0);
