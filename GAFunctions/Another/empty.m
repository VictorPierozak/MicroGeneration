function [z,ceq] = empty(x)
[row,~] = size(x);
z = zeros(row,1);
ceq = [];
end