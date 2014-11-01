%% Two state channel model with output Y given by (0/Unbound, 1/Bound) and
% input X given by (0/Low, 1/High).  
%
% Calculate the discrete time transition matrix, and bounds on the mutual
% information rates.  
%
% Background.  The mutual information rate is
%
% MI(X,Y) = H(X)+H(Y)-H(X,Y)
%
% We take X to be a Markov process with transition probabilities 
%
% X=0->X=1 w.pr. r (per time step)
% X=1->X=0 w.pr. s (per time step)
% 
% Therefore H(X) is a known function of r and s.
%
% The transition matrix for X is [r,1-r; s,1-s]. (Rows summing to one.)
%
% The steady state distribution for X is

px0=@(r,s)s./(r+s); % probability that X=0 (low concentration)
px1=@(r,s)r./(r+s); % probability that X=1 (high concentration)

% The binary entropy function is 

phi2=@(p)(p.*log2(1./p)+(1-p).*log2(1./(1-p))); % entropy is measured in bits

% The entropy rate of X is 

HX=@(r,s)px0(r,s).*phi2(r)+px1(r,s).*phi2(s); % standard entropy rate for 2-state Markov process

%% Check
rvec=linspace(0,1,201); % r-grid for plotting
svec=linspace(0,1,200); % s-grid for plotting
[rplot,splot]=meshgrid(rvec,svec);
HXplot=HX(rplot,splot);
figure
[c,h]=contour(rplot,splot,HXplot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)

%% The joint process. 
% (X-index shifted by 1 to match notation in our paper, which was chosen to
% match the notation in the Berger papers.  The joint Markov process is
% then (X(t+1),Y(t))->(X(t+2),Y(t+1)).
%
% For notation we use X=0(Low) or X=1(High) and Y=0(Unbound) or Y=1(Bound).
%
% The output Y has transition probability Y(t)=0->Y(t+1)=1 that depends on
% the value of X(t+1).  We assume the transition probability
% Y(t)=1->Y(t+1)=0 is independent of X(t+1).  In the general case, we have the following.
% If X(t+1)=0 then Pr(Y(t+1)=1 | Y(t)=0) = alo, and Pr(Y(t+1)=0 | Y(t)=0) = 1-alo.
% If X(t+1)=1 then Pr(Y(t+1)=1 | Y(t)=0) = ahi, and Pr(Y(t+1)=0 | Y(t)=0) = 1-ahi.
% For all X, Pr(Y(t+1)=0 | Y(t)=1) = b, and Pr(Y(t+1)=1 | Y(t)=1) = 1-b.
%
% With this definition, the combined process Z(t)=(X(t+1),Y(t)) is a Markov
% process on four states:
% Z1=(X=0,Y=0)
% Z2=(X=0,Y=1)
% Z3=(X=1,Y=0)
% Z4=(X=1,Y=1)
%
% The transition matrix for Z is given by T (rows sum to 1)
T=@(r,s,alo,ahi,b)...
    [(1-r).*(1-alo),  (1-r).*alo,     r.*(1-alo),       r.*alo;
    (1-r).*b,           (1-r).*(1-b),  r.*b,               r.*(1-b);
    s.*(1-ahi),        s.*ahi,          (1-s).*(1-ahi), (1-s).*ahi;
    s.*b,                s.*(1-b),       (1-s).*b,         (1-s).*(1-b)];

% For future reference:
aloplot=.1; % nominal "low" binding rate
ahiplot=.9; % nominal "high" binding rate
bplot=.5; % nominal unbinding rate

% We will need the components of T one by one later.

% First Row of T:
TLU2LU=@(r,s,alo,ahi,b)(1-r).*(1-alo);
TLU2LB=@(r,s,alo,ahi,b)(1-r).*alo;
TLU2HU=@(r,s,alo,ahi,b)r.*(1-alo);
TLU2HB=@(r,s,alo,ahi,b)r.*alo;

% Second Row of T:
TLB2LU=@(r,s,alo,ahi,b)(1-r).*b;
TLB2LB=@(r,s,alo,ahi,b)(1-r).*(1-b);
TLB2HU=@(r,s,alo,ahi,b)r.*b;
TLB2HB=@(r,s,alo,ahi,b)r.*(1-b);

% Third Row of T:
THU2LU=@(r,s,alo,ahi,b)s.*(1-ahi);
THU2LB=@(r,s,alo,ahi,b)s.*ahi;
THU2HU=@(r,s,alo,ahi,b)(1-s).*(1-ahi);
THU2HB=@(r,s,alo,ahi,b)(1-s).*ahi;

% Fourth Row of T:
THB2LU=@(r,s,alo,ahi,b)s.*b;
THB2LB=@(r,s,alo,ahi,b)s.*(1-b);
THB2HU=@(r,s,alo,ahi,b)(1-s).*b;
THB2HB=@(r,s,alo,ahi,b)(1-s).*(1-b);

% To check: alternate definition of T.  It works! (2013-03-15 2:00 pm EDT)
% Now comment it out and use the original only.
%T=@(r,s,alo,ahi,b)...
%    [TLU2LU(r,s,alo,ahi,b),TLU2LB(r,s,alo,ahi,b),TLU2HU(r,s,alo,ahi,b),TLU2HB(r,s,alo,ahi,b);...
%    TLB2LU(r,s,alo,ahi,b),TLB2LB(r,s,alo,ahi,b),TLB2HU(r,s,alo,ahi,b),TLB2HB(r,s,alo,ahi,b);...
%    THU2LU(r,s,alo,ahi,b),THU2LB(r,s,alo,ahi,b),THU2HU(r,s,alo,ahi,b),THU2HB(r,s,alo,ahi,b);...
%    THB2LU(r,s,alo,ahi,b),THB2LB(r,s,alo,ahi,b),THB2HU(r,s,alo,ahi,b),THB2HB(r,s,alo,ahi,b)];

%% The stationary distribution for Z is given by

abar=@(r,s,alo,ahi)(r*ahi + s*alo)./(r + s); % average of "a" over input ensemble
lambda=@(r,s)1-(r+s); % eigenvalue for convergence rate of 2 state Markov input process
% make sure r and s can be arrays, for plotting
% ZZ is the normalization constant for the stationary probability
% distribution.
ZZ=@(r,s,alo,ahi,b)...
   b*s.*(1-lambda(r,s)+(ahi+b)*lambda(r,s))+...
    s.*(lambda(r,s)*alo*(ahi+b)+(r+s).*abar(r,s,alo,ahi))+...
    b*r.*(1-lambda(r,s)+(alo+b)*lambda(r,s))+...
    r.*(lambda(r,s)*ahi*(alo+b)+(r+s).*abar(r,s,alo,ahi));

pxy00=@(r,s,alo,ahi,b)b*s.*(1-lambda(r,s)+(ahi+b)*lambda(r,s))./ZZ(r,s,alo,ahi,b);
pxy01=@(r,s,alo,ahi,b)s.*(lambda(r,s)*alo*(ahi+b)+(r+s).*abar(r,s,alo,ahi))./ZZ(r,s,alo,ahi,b);
pxy10=@(r,s,alo,ahi,b)b*r.*(1-lambda(r,s)+(alo+b)*lambda(r,s))./ZZ(r,s,alo,ahi,b);
pxy11=@(r,s,alo,ahi,b)r.*(lambda(r,s).*ahi*(alo+b)+(r+s).*abar(r,s,alo,ahi))./ZZ(r,s,alo,ahi,b);

%% test stationary distribution
commandwindow
disp('Yowza!')
for j=1:5
    r=rand;s=rand;alo=rand;ahi=rand;b=rand;
    pxy=[pxy00(r,s,alo,ahi,b),pxy01(r,s,alo,ahi,b),pxy10(r,s,alo,ahi,b),pxy11(r,s,alo,ahi,b)];
    disp(sum(pxy))
    disp([pxy;pxy*T(r,s,alo,ahi,b)]); % should be 1 and identical)
    pause(1)
end
clear pxy r s alo ahi b
% OK, at long last it works!  2013-03-15, 11:30 am EDT

%% Entropy rate of (X,Y) together 
% In the paper's notation, this is the entropy of the process
% Z(t)=(X(t+1),Y(t)).  

% We will need the entropy of a four-component distribution.  The
% quaternary entropy function is (in bits)

phi4=@(p1,p2,p3,p4)...
    (p1.*log2(1./p1)+p2.*log2(1./p2)+p3.*log2(1./p3)+p4.*log2(1./p4)); 

HXY=@(r,s,alo,ahi,b)...
    pxy00(r,s,alo,ahi,b).*...
    phi4(TLU2LU(r,s,alo,ahi,b),TLU2LB(r,s,alo,ahi,b),TLU2HU(r,s,alo,ahi,b),TLU2HB(r,s,alo,ahi,b))...
    +...
    pxy01(r,s,alo,ahi,b).*...
    phi4(TLB2LU(r,s,alo,ahi,b),TLB2LB(r,s,alo,ahi,b),TLB2HU(r,s,alo,ahi,b),TLB2HB(r,s,alo,ahi,b))...
    +...
    pxy10(r,s,alo,ahi,b).*...
    phi4(THU2LU(r,s,alo,ahi,b),THU2LB(r,s,alo,ahi,b),THU2HU(r,s,alo,ahi,b),THU2HB(r,s,alo,ahi,b))...
    +...
    pxy11(r,s,alo,ahi,b).*...
    phi4(THB2LU(r,s,alo,ahi,b),THB2LB(r,s,alo,ahi,b),THB2HU(r,s,alo,ahi,b),THB2HB(r,s,alo,ahi,b));

%   Entries of "T": each row is a probability distribution of destinations.        
%   [TLU2LU(r,s,alo,ahi,b),TLU2LB(r,s,alo,ahi,b),TLU2HU(r,s,alo,ahi,b),TLU2HB(r,s,alo,ahi,b);...
%   TLB2LU(r,s,alo,ahi,b),TLB2LB(r,s,alo,ahi,b),TLB2HU(r,s,alo,ahi,b),TLB2HB(r,s,alo,ahi,b);...
%   THU2LU(r,s,alo,ahi,b),THU2LB(r,s,alo,ahi,b),THU2HU(r,s,alo,ahi,b),THU2HB(r,s,alo,ahi,b);...
%   THB2LU(r,s,alo,ahi,b),THB2LB(r,s,alo,ahi,b),THB2HU(r,s,alo,ahi,b),THB2HB(r,s,alo,ahi,b)];

%% plot HXY

HXYplot=HXY(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,HXYplot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Joint Entropy Rate H(X,Y)','FontSize',20)

% Checks out.  Had to replace "*" with ".*" for r, s products to work
% properly.  PJT 2013-03-15 2:30 p.m.

%% The Mutual information is MI=H(X)+H(Y)-H(X,Y)
% Since H(Y)>=0, MI >=H(X)-H(X,Y) gives the simplest lower bound. Plot this
% lower bound "0"
MILB0=@(r,s,alo,ahi,b) HX(r,s)-HXY(r,s,alo,ahi,b);
MILB0plot=MILB0(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,MILB0plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Lower Bound 0: H(X)-H(X,Y)','FontSize',20)

%% Entropy of Y (single trial)
py0=@(r,s,alo,ahi,b)pxy00(r,s,alo,ahi,b)+pxy10(r,s,alo,ahi,b);
py1=@(r,s,alo,ahi,b)pxy01(r,s,alo,ahi,b)+pxy11(r,s,alo,ahi,b);
HY=@(r,s,alo,ahi,b)phi2(py0(r,s,alo,ahi,b));
HYcheck=@(r,s,alo,ahi,b)phi2(py1(r,s,alo,ahi,b)); % should be same as HY.

%Check construction
r=rand(1,5);s=rand(1,5);alo=rand;ahi=rand;b=rand;
commandwindow
disp([HY(r,s,alo,ahi,b);HYcheck(r,s,alo,ahi,b)])
clear r s alo ahi b

% plot HY
HYplot=HY(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,HYplot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
shg

% Looks OK (2014-04-06)

%% See Capacity2StateMarkov2.m, picking up from this point, with more compact notation!
% 2014-04-06

%% Calculate HYgX=H[Y_{t+1}|X_{t+1}]
% conditional probabilities
% Notation:
% Pr[Y_{t+1}=U|X_{t+1}=L], for example, is represented by this variable:
% pytp1eugxtp1el, and so on. Cumbersome, but unambiguous!
% First need these joint probabilities (see Table 4 of
% Capacity2StateMarkov) of Pr[X_{t+1},Y_{t+1}]
pxtp1elaytp1eu=@(r,s,alo,ahi,b)(1-alo)*pxy00(r,s,alo,ahi,b)+b*pxy01(r,s,alo,ahi,b);
pxtp1elaytp1eb=@(r,s,alo,ahi,b)alo*pxy00(r,s,alo,ahi,b)+(1-b)*pxy01(r,s,alo,ahi,b);
pxtp1ehaytp1eu=@(r,s,alo,ahi,b)(1-ahi)*pxy10(r,s,alo,ahi,b)+b*pxy11(r,s,alo,ahi,b);
pxtp1ehaytp1eb=@(r,s,alo,ahi,b)ahi*pxy10(r,s,alo,ahi,b)+(1-b)*pxy11(r,s,alo,ahi,b);
% Conditional probabilities from Bayes' formula
pytp1eugxtp1el=@(r,s,alo,ahi,b)pxtp1elaytp1eu(r,s,alo,ahi,b)./...
    (pxtp1elaytp1eu(r,s,alo,ahi,b)+pxtp1elaytp1eb(r,s,alo,ahi,b));
pytp1eugxtp1eh=@(r,s,alo,ahi,b)pxtp1ehaytp1eu(r,s,alo,ahi,b)./...
    (pxtp1ehaytp1eu(r,s,alo,ahi,b)+pxtp1ehaytp1eb(r,s,alo,ahi,b));
% Now the conditional entropy
HYgX=@(r,s,alo,ahi,b)...
    phi2(pytp1eugxtp1el(r,s,alo,ahi,b)).*px0(r,s)+...
    phi2(pytp1eugxtp1eh(r,s,alo,ahi,b)).*px1(r,s);

%% Calculate HYgYX=H[Y_{t+1}|X_{t+1},Y_t]
% from p. 22 of Capacity2StateMarkov.pdf
HYgYX=@(r,s,alo,ahi,b)...
    phi2(alo)*pxy00(r,s,alo,ahi,b)+...
    phi2(ahi)*pxy10(r,s,alo,ahi,b)+...
    phi2(b)*(pxy01(r,s,alo,ahi,b)+pxy11(r,s,alo,ahi,b));

HYgYXplot=HYgYX(rplot,splot,aloplot,ahiplot,bplot);

%% Calculate HYgYYX=H[Y_{t+2}|X_{t+1},Y_t,Y_{t+1}]
% which is the same thing as H[Y_{t+2}|X_{t+1},Y_{t+1}} since Y_{t+2} is
% independent of Y_t, given Y_{t+1} and X_{t+1}.
HYgYYX=@(r,s,alo,ahi,b)...
    (pxy00(r,s,alo,ahi,b)*(1-alo)+pxy01(r,s,alo,ahi,b)*b).*phi2((1-r)*alo+r*ahi)+...
    (pxy00(r,s,alo,ahi,b)*alo+pxy01(r,s,alo,ahi,b)*(1-b)).*phi2(b)+...
    (pxy10(r,s,alo,ahi,b)*(1-ahi)+pxy11(r,s,alo,ahi,b)*b).*phi2(s*alo+(1-s)*ahi)+...
    (pxy10(r,s,alo,ahi,b)*ahi+pxy11(r,s,alo,ahi,b)*(1-b)).*phi2(1-b);

%% MI Approximations

MIY=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HY(r,s,alo,ahi,b);
MIYgX=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYgX(r,s,alo,ahi,b);
MIYgYX=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYgYX(r,s,alo,ahi,b);
MIYgYYX=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYgYYX(r,s,alo,ahi,b);

%% Plot MI Approximations
MIYplot=MIY(rplot,splot,aloplot,ahiplot,bplot);
MIYgXplot=MIYgX(rplot,splot,aloplot,ahiplot,bplot);
MIYgYXplot=MIYgYX(rplot,splot,aloplot,ahiplot,bplot);
MIYgYYXplot=MIYgYYX(rplot,splot,aloplot,ahiplot,bplot);

figure
surf(rplot,splot,MIYplot); shading flat
text(.5,.999,MIY(.5,1,aloplot,ahiplot,bplot),'MIY','FontSize',20)
hold on
surf(rplot,splot,MIYgYYXplot); shading flat
text(.5,.999,MIYgYYX(.5,1,aloplot,ahiplot,bplot),'MIYgYYX')
surf(rplot,splot,MIYgYXplot); shading flat
text(.5,.999,MIYgYX(.5,1,aloplot,ahiplot,bplot),'MIYgYX')
surf(rplot,splot,MIYgXplot); shading flat
text(.5,.999,MIYgX(.5,1,aloplot,ahiplot,bplot),'MIYgX')

set(gca,'FontSize',20)
xlabel('r','FontSize',20),ylabel('s','FontSize',20)
view([160,20])
rotate3d on
shg
%print -dpdf MIapprox_layers

%% Result from the 2013 ISIT submission.
% Andrew Eckford (using results from Chen and Berger, and Berger and Ying)
% showed that the capacity is 

C_Eckford=@(r,s,alo,ahi,b) ...
    (phi2(ahi*px1(r,s)+alo*px0(r,s))...
    -px1(r,s).*phi2(ahi)...
    -px0(r,s).*phi2(alo))./...
    (1+(ahi*px1(r,s)+alo*px0(r,s))./b);

%This result holds when the input is IID, i.e. when r+s=1, only.  The
%commented out material below plots the expression in the entire (r,s)
%plane.  Instead, let's compare C_Eckford with the entropy rate for the
%input (which should be an upper bound) along the line s=1-r.

figure
rplot1=0:.01:1;
splot1=1-rplot1;
plot(rplot1,C_Eckford(rplot1,splot1,aloplot,ahiplot,bplot),'LineWidth',3)
hold on
plot(rplot1,HX(rplot1,splot1),'g-','LineWidth',3) 
set(gca,'FontSize',20)
xlabel('Up Rate r','FontSize',20)
ylabel('Information Rate','FontSize',20)
legend('MI rate','Input rate')
title('MI rate for \alpha_-=0.1, \alpha_+=0.9, \beta=0.5','FontSize',20)
set(gcf,'PaperOrientation','landscape')
print -dpdf C_ISIT2013.pdf
shg

% figure
% 
% subplot(2,2,1);Cap_plot=C_Eckford(rplot,splot,aloplot,ahiplot,bplot);
% [c2,h2]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.1,\alpha_+=0.9,\beta=0.5','FontSize',20)
% clabel(c2,h2);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% 
% subplot(2,2,2);Cap_plot=C_Eckford(rplot,splot,.1,.9,.9);
% [c3,h3]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.1,\alpha_+=0.9,\beta=0.9','FontSize',20)
% clabel(c3,h3);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% 
% subplot(2,2,3);Cap_plot=C_Eckford(rplot,splot,.01,.99,.5);
% [c4,h4]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.01,\alpha_+=0.99,\beta=0.5','FontSize',20)
% clabel(c4,h4);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% 
% subplot(2,2,4);Cap_plot=C_Eckford(rplot,splot,.01,.99,.99);
% [c5,h5]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.01,\alpha_+=0.99,\beta=0.99','FontSize',20)
% clabel(c5,h5);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% 
% set(gcf,'PaperOrientation','landscape')
% print -dpdf C_ISIT2013.pdf
% 
% %% same thing, zoomed in
% 
% figure
% 
% subplot(2,2,1);Cap_plot=C_Eckford(rplot,splot,aloplot,ahiplot,bplot);
% [c2,h2]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.1,\alpha_+=0.9,\beta=0.5','FontSize',20)
% clabel(c2,h2);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% axis([0 .1 0 .1])
% 
% subplot(2,2,2);Cap_plot=C_Eckford(rplot,splot,.1,.9,.9);
% [c3,h3]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.1,\alpha_+=0.9,\beta=0.9','FontSize',20)
% clabel(c3,h3);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% axis([0 .1 0 .1])
% 
% subplot(2,2,3);Cap_plot=C_Eckford(rplot,splot,.01,.99,.5);
% [c4,h4]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.01,\alpha_+=0.99,\beta=0.5','FontSize',20)
% clabel(c4,h4);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% axis([0 .1 0 .1])
% 
% subplot(2,2,4);Cap_plot=C_Eckford(rplot,splot,.01,.99,.99);
% [c5,h5]=contour(rplot,splot,Cap_plot);title('\alpha_-=0.01,\alpha_+=0.99,\beta=0.99','FontSize',20)
% clabel(c5,h5);axis equal,axis tight,axis([0 1 0 1]);colorbar;set(gca,'FontSize',20);shg
% xlabel('r','FontSize',20),ylabel('s','FontSize',20)
% axis([0 .1 0 .1])

%% 