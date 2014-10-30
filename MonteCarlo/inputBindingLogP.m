function y = inputBindingLogP(x,b,P,PB,py,px)

%sx = size(P);
%numStates = sx(1);

if (nargin == 5)
    foo = P^1000;
    px = foo(1,:);
end

result = log(px(x(1))*py(b(1)));

for c = 2:length(x)
    tempP = P(x(c-1),x(c));
    tempP = tempP * PB(b(c-1),b(c),x(c-1));
    result = result + log(tempP);
end

y = result;