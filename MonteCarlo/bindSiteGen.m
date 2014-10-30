function y = bindSiteGen(x,PB,py)

% x = 1 --> low, x = 2 --> high
% py is the prior on binding site probability

result = zeros(1,length(x));

result(1) = randNonUniform(1,1,py);

for c = 2:length(x)
    result(c) = randNonUniform(1,1,PB(result(c-1),:,x(c-1)));
end


y = result;