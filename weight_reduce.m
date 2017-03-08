%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 10000;

% A = HW
dictionary = load('dictionary/dict50.mat');
dictionary = dictionary.dictionary;
dictionary = real(dictionary);

H = load('dictionary/H50.mat');
H = H.H;

W = load('dictionary/W50.mat');
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

WEIGHT = [];

tic
[YFIT0,R,COEFF,IOPT] = wmpalg('BMP', signal, W');
%disp(IOPT);
nIOPT = unique(IOPT);
for n=1:length(nIOPT)
	idx = IOPT==nIOPT(n);
	w = idx*COEFF;
	WEIGHT = [WEIGHT w];
end
disp(WEIGHT);
Ht = H';
%new_dict = dictionary(find(abs(sum(Ht(nIOPT',:)))>0.8),:);
tmp = Ht(nIOPT',:);
Max = max(abs(WEIGHT*tmp));
new_dict = dictionary(find(abs(WEIGHT*tmp)>0.3*Max),:);
disp(size(new_dict));
[YFIT1,R,COEFF,IOPT0] = wmpalg('BMP', signal, new_dict');
toc

audiowrite('new.wav', YFIT0, Fs);
audiowrite('refine.wav', YFIT1, Fs);
err1 = sum(abs(YFIT0-signal'));
err2 = sum(abs(YFIT1-signal'));

disp(max(abs(WEIGHT*tmp)));

disp(err0);
disp(err1);
disp(err2);
disp((err1-err0)/err0);
disp((err2-err0)/err0);
