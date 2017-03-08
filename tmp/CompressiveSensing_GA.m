[x, fval] = ga(@fitness, 1000);
siglen = 200;
seq = 1:siglen;
signal = sin(seq).*(cos(log(seq)));
tic
%[x,fval] = gamultiobj(@fitness, 1000);
toc
disp(fval);
disp(sum(abs(x-signal')));