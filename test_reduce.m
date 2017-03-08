%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 10000;

% A = HW => A' = W'H'
dictionary = load('dictionary/dict60.mat');
dictionary = dictionary.dictionary;
dictionary = real(dictionary);

H = load('dictionary/H60_100.mat');
H = H.H;

W = load('dictionary/W60_100.mat');
W = W.W;

[signal, Fs] = audioread('whale.wav');
signal = signal';


disp(size(dictionary));
tic
[YFIT,R,COEFF, iopt] = wmpalg('BMP', signal, dictionary');
toc

audiowrite('original.wav', YFIT, Fs);
err0 = sum(abs(YFIT-signal'));

disp(size(H));
disp(size(W));

tic
[YFIT0,R,COEFF,IOPT] = wmpalg('BMP', signal, W');
IOPT = unique(IOPT);
Ht = H';

tmp = Ht(IOPT,:);
Max = max(abs(sum(Ht(IOPT',:))));
new_dict = dictionary(find(abs(sum(Ht(IOPT,:)))>0.78*Max),:);
[YFIT1,R,COEFF,IOPT0] = wmpalg('BMP', signal, new_dict');
toc

audiowrite('new.wav', YFIT0, Fs);
audiowrite('refine.wav', YFIT1, Fs);
err1 = sum(abs(YFIT0-signal'));
err2 = sum(abs(YFIT1-signal'));

disp(size(new_dict));
disp(size(sum(Ht(IOPT',:))));
disp(max(abs(sum(Ht(IOPT',:)))));

disp(err0);
disp(err1);
disp(err2);
disp((err1-err0)/err0);
disp((err2-err0)/err0);
