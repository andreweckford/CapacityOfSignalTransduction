clear
load result
load result2
load result3
load result4
z = z1+z2+z3+z4;

map = colormap;
map(:,3) = (128:-0.5:96.5)./128;
map(:,1) = (65:128)./128;
map(:,2) = map(:,1);
colormap(map);
[c,h] = contourf(ss,rr,z); 
axis square;
xlabel('s (high-to-low transition probability)')
ylabel('r (low-to-high transition probability)')
title('\alpha_H = 0.5, \alpha_L = 0.1, \beta = 0.2')
clabel(c,h,'manual');
