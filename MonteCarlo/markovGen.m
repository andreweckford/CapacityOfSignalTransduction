function y = markovGen(len,P,ps)

[rows cols] = size(P);
S = 1:cols; % state alphabet
result = zeros(1,len);

if (nargin == 2)
    % a smarter way to do this, find the eigenvectors
    foo = P^1000;
    ps = foo(1,:);
end

result(1) = randNonUniform(1,1,ps);
for c=2:len
    result(c) = randNonUniform(1,1,P(result(c-1),:));
end

y = result;