function y = mutualInfoMC(P,PB,py,len,iter)

result = 0;

for c = 1:iter
    
    x = markovGen(len,P);
    b = bindSiteGen(x,PB,py);
    
    foo = inputBindingLogP(x,b,P,PB,[1 0]) ...
        - inputLogP(x,P) ...
        - bindingLogP(b,P,PB,[1 0]);
    
    result = result + foo;
    
    %if (foo < 0) warning('Houston, we have a problem'), end
    
end

y = result/len/iter/log(2); 