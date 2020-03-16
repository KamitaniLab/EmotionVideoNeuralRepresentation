function y = entropy0(x)
% entropy -- calcurate enctropy
% y = entropy0(x)
%
% [Input]
%  x: probability vector or matrix
%
% [Output]
%  y: entropy
%
% [Usage]
%  entropy0([1/2 1/4 1/8 1/16 1/64 1/64 1/64 1/64]) ; see PRML chap1.6
%
% [note]
% 
% written by: Tomoyasu Horikawa  horikawa-t@atr.jp 2019/08/24
%
y = 0;
x(x(:) == 0) = [];
x = x./sum(x);
for i = 1:length(x)
    e = x(i) * -log2(x(i));
    if isnan(e)
        afaf
    else
        y = y + e;
    end
end
