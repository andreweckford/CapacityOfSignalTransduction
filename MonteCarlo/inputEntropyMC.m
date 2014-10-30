function y = inputEntropyMC(P,len,iter)

result = 0;

for c = 1:iter
    
    x = markovGen(len,P);
    
    foo = inputLogP(x,P);
    
    result = result - foo;
        
end

y = result/len/iter/log(2); 