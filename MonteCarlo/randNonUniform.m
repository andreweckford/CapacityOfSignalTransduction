function y = randNonUniform(rows,cols,p);

result = ones(rows,cols);
cp = cumsum(p);
r = rand(rows,cols);

for c=1:(length(p)-1)
    result = result + (r > cp(c));
end

y = result;