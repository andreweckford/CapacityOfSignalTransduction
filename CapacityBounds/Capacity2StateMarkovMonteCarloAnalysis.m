%% load data from the run
clear all, close all
cd ~/projects/bioinfo/code/pjthomas/matlab/capacity_binding
load Capacity2StateMarkovMonteCarloData2
nstep=1e6;

%% Empirical probability that X=1
px=sum(x,3)/nstep;

%% Empirical probability that Y=1
py=sum(y,3)/nstep;

%% Empirical probability that (X,Y)=(0,0) or (0,1) or (1,0) or (1,1) (in
% that order)
xy=2*x+y; % binary word such that (0,0)=>0, (0,1)=>1, (1,0)=>2, (1,1)=>3.
for j=1:4,
    pxy(:,:,j)=sum((xy==j-1),3)/nstep;
end

%% Empirical probability Pr(2*Y_{t+1}+Y_t)=0,1,2, or 3, is 
for j=1:4
    pyy(:,:,j)=sum(1+2*y(:,:,2:end-1)+y(:,:,1:end-2)==j,3)/(nstep-1);
end

%% Empirical probability Pr(Y_{t+2}=1|2*Y_{t+1}+Y_t=j), j \in {0,1,2,3} is
for j=1:4
    pygyy(:,:,j)=sum(y(:,:,3:end).*(1+2*y(:,:,2:end-1)+y(:,:,1:end-2)==j),3)./...
        (pyy(:,:,j)*(nstep-1));
end

%% Empirical probability Pr(4*X_t+2*Y_{t+1}+Y_t)
for j=1:8
    pxyy(:,:,j)=sum(1+4*x(:,:,1:end-2)+2*y(:,:,2:end-1)+y(:,:,1:end-2)==j,3)/(nstep-2);
end

%% Empirical probability Pr(Y_{t+2}=1|4*X_t+2*Y_{t+1}+Y_t), j \in {0,1,...,7} is
for j=1:8
    pygxyy(:,:,j)=sum(y(:,:,3:end).*(1+4*x(:,:,1:end-2)+2*y(:,:,2:end-1)+y(:,:,1:end-2)==j),3)./...
        (pxyy(:,:,j)*(nstep-2));
end

%% 2-fold entropy function
phi2=@(p)p.*log(1./p)+(1-p).*log(1./(1-p));

%% n-fold entropy function (assumes p is nr x ns x nconditions
phin=@(p)sum(p.*log(1./p),3);

%% Plot the entropy of X (zero step)
figure
pcolor(r,s,phi2(px))
shading flat, colormap gray,axis equal,axis tight,title('H[X]'),shg

%% Plot the entropy of Y (zero step)
figure
pcolor(r,s,phi2(py))
shading flat, colormap gray,axis equal,axis tight,title('H[Y]'),shg

%% Plot the entropy of (X,Y) (joint, zero step)
figure
pcolor(r,s,phin(pxy))
shading flat, colormap gray,axis equal,axis tight,title('H[Y]'),shg

%% Plot the H(Y) (zero step) MI approximation
figure
pcolor(r,s,phi2(px)+phi2(py)-phin(pxy))
shading flat, colormap gray,axis equal,axis tight,title('MI approx 0'),shg

%% Empirical entropy H[Y_{t+2}|Y_{t+1},Y_t]
HYgYY=zeros(ns,nr);
for j=1:4
    HYgYY=HYgYY+phi2(pygyy(:,:,j)).*pyy(:,:,j);
end
HYgYY=squeeze(HYgYY);
figure,surf(r,s,HYgYY),shading flat, shg

