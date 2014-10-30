function y = inputLogP(s, P, ps)

if (nargin == 2)
    foo = P^1000;
    ps = foo(1,:);
end

result = log(ps(s(1)));
for c = 2:length(s)
   result = result + log(P(s(c-1),s(c)));
end

y = result;