function [x,y,pxHyU,pyU] = genXY(PB,Pxy,len,px,py)

% Pxy: matrix of pr(x given y)
% first row: y = UU(11); second row: y = UB(12);
% third row: Y = BU(21); fourth row: y = BB(22) (rightmost is most recent)
% first column: x = L(1); second column: x = H(2)

x = zeros(1,len);
y = zeros(1,len);

pxHyU = 0;
pyU = 0;
y(1) = randNonUniform(1,1,py);
y(2) = randNonUniform(1,1,py);
x(1) = randNonUniform(1,1,px);
x(2) = randNonUniform(1,1,px);

for c = 3:len
    
    % feedback x
    yIndex = y(c-1)-1 + 2*(y(c-2)-1) + 1;
    x(c) = randNonUniform(1,1,Pxy(yIndex,:));
    y(c) = randNonUniform(1,1,PB(y(c-1),:,x(c)));
    if ((x(c) == 2)&&(y(c-1) == 1))
        pxHyU = pxHyU + 1;
    end
    if (y(c) == 1)
        pyU = pyU + 1;
    end
        
end

pxHyU = pxHyU / pyU;
pyU = pyU / len;