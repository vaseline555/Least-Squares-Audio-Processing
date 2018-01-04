function [x] = back_substitution(R, QTb)
%%% From Gram-Schmidt algorithm, get QR-factorization, using Q and R to
%%% solve Ax = b <=> QRx = b <=> Rx = (Q^T)*b
%%% For Rx = (Q^T)*b, use this function.
n = length(QTb);
x(n, 1) = QTb(n) / R(n,n);
for i = n - 1 : -1 : 1
    x(i, 1) = (QTb(i) - R(i, i + 1 : n) * x(i + 1 : n, 1)) ./ R(i, i);
end