function y = hCondXbrute(len,trials,P,PB,px,py)

% Calculates the entropy of x_i given x^{i-1},y^{i-1}
% Calculation is done "the hard way"

result = 0;

for c = 1:trials
    x = markovGen(len,P,px);
    y = bindSiteGen(x,PB,py);
    y(len) = 1;
    num = exp(inputBindingLogP(x,y,P,PB,py,px));
    y(len) = 2;
    num = num + exp(inputBindingLogP(x,y,P,PB,py,px));
    den = exp(inputBindingLogP(x(1:len-1),y(1:len-1),P,PB,py,px));
    
    result = result - log2(num/den);
end

y = result;