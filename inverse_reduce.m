%%% it is used to test matchin pursuit with NMF reduction%%%
% parameters: t0, f0, sigma
siglen = 10000;

% A = HW, A' = W'H'
dictionary = load('unique_dict/dict1000_97.mat');
dictionary = dictionary.dictionary;
dictionary = real(dictionary);

H = load('unique_dict/H1000_97_300.mat');
H = H.H;

W = load('unique_dict/W1000_97_300.mat');
W = W.W;

Sd = size(dictionary);

[signal, Fs] = audioread('whale.wav');
signal = signal';

disp(size(dictionary));
tic
[YFIT,R,COEFF, iopt] = wmpalg('BMP', signal, dictionary');
toc

iopt = unique(iopt);
audiowrite('original.wav', YFIT, Fs);
err0 = sum(abs(YFIT-signal'));

disp(size(H));
disp(size(W));

WEIGHT = [];
Ht = H';
invH = pinv(Ht);
allweight = zeros(1,size(Ht,1));
oneweight = zeros(1,size(Ht,1));

tic
[YFIT0,R,COEFF,IOPT] = wmpalg('BMP', signal, W');
nIOPT = unique(IOPT);
for n=1:length(nIOPT)
	idx = IOPT==nIOPT(n);
	w = idx*COEFF;
	WEIGHT = [WEIGHT w];
end
fprintf('weight of matching pursuit algorithm.\n');
disp(WEIGHT);
allweight(nIOPT) = WEIGHT;
oneweight(nIOPT) = 1;
%new_dict = dictionary(find(abs(sum(Ht(nIOPT',:)))>0.8),:);
%disp(size(invH));
%disp(size(allweight));
tmp = invH*allweight';
Max = max(abs(tmp));
new_dict = dictionary(find(abs(tmp)>0.01*Max),:);

fprintf('check if true atoms are in the new dictionary:\n');
fprintf('Maximum value: %4.8f\n', Max);
fprintf('ratio to MAX for useful atoms\n:');
disp(abs(tmp(iopt))'/Max);
%disp(size(tmp(find((abs(tmp)/Max)>0.05))));
fprintf('mean ratio for all atoms: %3.4f\n', mean(abs(tmp)/Max));


%fprintf('abs value weight gain\n');
%tmp2 = abs(invH)*oneweight';
%Max = max(abs(tmp2));
%disp(abs(tmp2(iopt))'/Max);
%fprintf('number of atoms with ratio larger than %3.4f\n', 0.05);
%disp(size(tmp2(find((abs(tmp2)/Max)>0.05))));
%%%%%%%%%%%%%%%%%%%%%%%%
%new_dict = dictionary(find(abs(tmp)>0.00*Max),:);
%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('mean ratio:');
%disp(mean(abs(tmp2))/Max);

d = size(new_dict, 1);
fprintf('the size of dictionary for proposed method:\n');
disp(size(new_dict));
[YFIT1,R,COEFF,IOPT0] = wmpalg('BMP', signal, new_dict');
toc


%audiowrite('new.wav', YFIT0, Fs);
%audiowrite('refine.wav', YFIT1, Fs);
err1 = sum(abs(YFIT0-signal'));
err2 = sum(abs(YFIT1-signal'));

tic
idx = randperm(Sd(1),d);
rdict = dictionary(idx,:);
fprintf('the size of dictionary for random sampling method:\n');
disp(size(rdict));
[YFIT0, R, COEFF, iopt1] = wmpalg('BMP', signal, rdict');
toc

err3 = sum(abs(YFIT0-signal'));

disp(err0);
%disp(err1);
r2 = (err2-err0)/err0;
fprintf('error and ratio of selected atoms by NMF: %3.4f, %3.4f\n', err2, r2)
r3 = (err3-err0)/err0;
fprintf('error and ratio of random selected atoms: %3.4f, %3.4f\n', err3, r3);
