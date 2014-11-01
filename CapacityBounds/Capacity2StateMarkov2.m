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

% Transition probability matrix in a different form:
% Pr[X3=x3,Y2=y2|X2=x2,Y1=y1], i.e. T(x2,y1->x3,y2), written as a function of x2,y1,x3,y2:
TT=@(r,s,alo,ahi,b,x2,y1,x3,y2)...
    ...% First Row of T:
    (1-x2).*(1-y1).*(1-x3).*(1-y2).*TLU2LU(r,s,alo,ahi,b)+...
    (1-x2).*(1-y1).*(1-x3).*y2.*TLU2LB(r,s,alo,ahi,b)+...
    (1-x2).*(1-y1).*x3.*(1-y2).*TLU2HU(r,s,alo,ahi,b)+...
    (1-x2).*(1-y1).*x3.*y2.*TLU2HB(r,s,alo,ahi,b)+...
    ...% Second Row of T:
    (1-x2).*y1.*(1-x3).*(1-y2).*TLB2LU(r,s,alo,ahi,b)+...
    (1-x2).*y1.*(1-x3).*y2.*TLB2LB(r,s,alo,ahi,b)+...
    (1-x2).*y1.*x3.*(1-y2).*TLB2HU(r,s,alo,ahi,b)+...
    (1-x2).*y1.*x3.*y2.*TLB2HB(r,s,alo,ahi,b)+...
    ...% Third Row of T:
    x2.*(1-y1).*(1-x3).*(1-y2).*THU2LU(r,s,alo,ahi,b)+...
    x2.*(1-y1).*(1-x3).*y2.*THU2LB(r,s,alo,ahi,b)+...
    x2.*(1-y1).*x3.*(1-y2).*THU2HU(r,s,alo,ahi,b)+...
    x2.*(1-y1).*x3.*y2.*THU2HB(r,s,alo,ahi,b)+...
    ...% Fourth Row of T:
    x2.*y1.*(1-x3).*(1-y2).*THB2LU(r,s,alo,ahi,b)+...
    x2.*y1.*(1-x3).*y2.*THB2LB(r,s,alo,ahi,b)+...
    x2.*y1.*x3.*(1-y2).*THB2HU(r,s,alo,ahi,b)+...
    x2.*y1.*x3.*y2.*THB2HB(r,s,alo,ahi,b);

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
disp('Checking stationary distribution code is consistent...')
for j=1:5
    r=rand;s=rand;alo=rand;ahi=rand;b=rand;
    pxy=[pxy00(r,s,alo,ahi,b),pxy01(r,s,alo,ahi,b),pxy10(r,s,alo,ahi,b),pxy11(r,s,alo,ahi,b)];
    disp(sum(pxy))
    disp([pxy;pxy*T(r,s,alo,ahi,b)]); % should be 1 and identical)
    pause(1)
end
clear pxy r s alo ahi b
disp('OK, at long last it works!  2013-03-15, 11:30 am EDT')

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

% Comment: this quantity turns out to be negative.  Since MI>=0, this is
% not much of a lower bound!

%% Entropy of Y (single trial) Gives the simplest upper bound on H[Y]
py0=@(r,s,alo,ahi,b)pxy00(r,s,alo,ahi,b)+pxy10(r,s,alo,ahi,b);
py1=@(r,s,alo,ahi,b)pxy01(r,s,alo,ahi,b)+pxy11(r,s,alo,ahi,b);
HYUB0=@(r,s,alo,ahi,b)phi2(py0(r,s,alo,ahi,b));
HYcheck=@(r,s,alo,ahi,b)phi2(py1(r,s,alo,ahi,b)); % should be same as HY.

%Check construction
disp('Quick check that I constructed the entropy H[Y] (one step) correctly:')
r=rand(1,5);s=rand(1,5);alo=rand;ahi=rand;b=rand;
commandwindow
disp([HYUB0(r,s,alo,ahi,b);HYcheck(r,s,alo,ahi,b)])
clear r s alo ahi b
disp('OK, appears consistent.  Moving on...')

% plot HY
HYplot=HYUB0(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,HYplot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
shg

% Looks OK (2014-04-06)

%% Simplest upper bound is MI<=H(X)-H(X,Y)+H(Y,snapshot)
MIUB0=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYUB0(r,s,alo,ahi,b);
MIUB0plot=MIUB0(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,MIUB0plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Upper Bound 0: H(X)-H(X,Y)+H(Y1)','FontSize',20)
shg

% Comment: this quantity shows the famous "ridge" that apparently stays
% positive all the way to r=s=0, where we know MI must go to zero.  

%% Next simplest Lower Bound is based on H[Y2|X2,Y1]
HYLB1=@(r,s,alo,ahi,b) ...
    pxy00(r,s,alo,ahi,b).*phi2(TLU2LB(r,s,alo,ahi,b)+TLU2HB(r,s,alo,ahi,b))+...
    pxy01(r,s,alo,ahi,b).*phi2(TLB2LB(r,s,alo,ahi,b)+TLB2HB(r,s,alo,ahi,b))+...
    pxy10(r,s,alo,ahi,b).*phi2(THU2LB(r,s,alo,ahi,b)+THU2HB(r,s,alo,ahi,b))+...
    pxy11(r,s,alo,ahi,b).*phi2(THB2LB(r,s,alo,ahi,b)+THB2HB(r,s,alo,ahi,b));
MILB1=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYLB1(r,s,alo,ahi,b);
MILB1plot=MILB1(rplot,splot,aloplot,ahiplot,bplot);
% Skip this plot -- MI is zero to within machine precision, but takes a
% long time to plot on account of the many interwoven level curves.
%figure
%[c,h]=contour(rplot,splot,MILB1plot);
%clabel(c,h);
%axis equal,axis tight,axis([0 1 0 1])
%colorbar
%set(gca,'FontSize',20)
%title('Lower Bound 1: H(X)-H(X,Y)+H(Y2|X2,Y1)','FontSize',20)
%shg
    
% Comment: this quantity apparently evaluates to zero, within machine
% precision!  So that is the same as saying MI>=0.

%% Next simplest Upper Bound is based on H[Y2|Y1]

% We need Pr[X3=x3,Y2=y2|Y1=y1] and Pr[x2=x2|Y1=y1].  

% To convert from (x,y) indexing to i=[1,2,3,4] indexing, note that
% [0:3=(2x+y)] has a one in exactly location 1,2,3 or 4, depending on
% whether xy=00, xy=01, xy=10, or xy=11, respectively.

% Pr[X2=x2|Y1=y1]
px2gy1=@(r,s,alo,ahi,b,x2,y1)...
    (1-x2).*(1-y1).*(pxy00(r,s,alo,ahi,b)./(py0(r,s,alo,ahi,b)))+... % case x2=0,y1=0
    (1-x2).*y1    .*(pxy01(r,s,alo,ahi,b)./(py1(r,s,alo,ahi,b)))+... % case x2=0,y1=1
    x2.*(1-y1)    .*(pxy10(r,s,alo,ahi,b)./(py0(r,s,alo,ahi,b)))+... % case x2=1,y1=0
    x2.*y1        .*(pxy11(r,s,alo,ahi,b)./(py1(r,s,alo,ahi,b)));    % case x2=1,y2=1

% Pr[X3=x3,Y2=y2|Y1=y1]
px3y2gy1=@(r,s,alo,ahi,b,x3,y2,y1)...
    TT(r,s,alo,ahi,b,0,y1,x3,y2)... % Pulls out the (L,y1)->(x3,y2) element of T, i.e. Pr[X3=x3,Y2=y2|X2=L,Y1=y1]
    .*px2gy1(r,s,alo,ahi,b,0,y1)+... % Gives Pr[X=L|Y1=y1]
    TT(r,s,alo,ahi,b,1,y1,x3,y2)... % Pulls out the (H,y1)->(x3,y2) element of T, i.e. Pr[X3=x3,Y2=y2|X2=H,Y1=y1]
    .*px2gy1(r,s,alo,ahi,b,1,y1);     % Gives Pr[X=H|Y1=y1]

HYUB1=@(r,s,alo,ahi,b) ...
    py0(r,s,alo,ahi,b).*phi2(px3y2gy1(r,s,alo,ahi,b,0,1,0)+px3y2gy1(r,s,alo,ahi,b,1,1,0))+...
    py1(r,s,alo,ahi,b).*phi2(px3y2gy1(r,s,alo,ahi,b,0,1,1)+px3y2gy1(r,s,alo,ahi,b,1,1,1));
MIUB1=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYUB1(r,s,alo,ahi,b);
MIUB1plot=min(MIUB1(rplot,splot,aloplot,ahiplot,bplot),HX(rplot,splot));
figure
[c,h]=contour(rplot,splot,MIUB1plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Upper Bound 1: min(H(X)-H(X,Y)+H(Y2|Y1),H(X))    ','FontSize',20)
shg
print -dpdf Cap2StateMarkov2_UB1.pdf

%% Compare the bounds we now have

bound_diff_1=@(r,s,alo,ahi,b)min(HX(r,s),MIUB1(r,s,alo,ahi,b))-MILB1(r,s,alo,ahi,b);
bd1plot=bound_diff_1(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,bd1plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Difference between UB1 and LB1  ','FontSize',20)
shg

%% Next Simplest Upper Bound is H[Y3|Y2,Y1]
% for this bound we need
% Pr[X3=x3,X2=x2|Y2=y2,Y1=y1]=px3x2gy2y1(x3,x2,y2,y1)
% and
% Pr[X4=x4,Y3=B|X3=x3,Y2=y2,X2=x2,Y1=y1]=px4y3gx3y2y1(x4,y3,x3,y2,y1), at
% Y3=B
% and 
% Pr[X3=x3,Y2=y2,X2=x2,Y1=y1]=px3y2x2y1(x3,y2,x2,y1).
% Then, H[Y3|Y1,Y2] is the sum over y1 and y2 of (
% sum over x2, x3 of (
% px3y2x2y1) )
% times 
% phi2(sum over x4 of (sum over x2,x3 of (prx4y3gx3y2x2y1*prx3x2gy2y1)))

py1x2y2x3=@(r,s,alo,ahi,b,y1,x2,y2,x3)(...
    (1-x2).*(1-y1).*pxy00(r,s,alo,ahi,b)+...
    (1-x2).*y1.*pxy01(r,s,alo,ahi,b)+...
    x2.*(1-y1).*pxy10(r,s,alo,ahi,b)+...
    x2.*y1.*pxy11(r,s,alo,ahi,b)).*...
    TT(r,s,alo,ahi,b,x2,y1,x3,y2);

py1y2=@(r,s,alo,ahi,b,y1,y2)(... % sum over x2,x3=0,1
    py1x2y2x3(r,s,alo,ahi,b,y1,0,y2,0)+...
    py1x2y2x3(r,s,alo,ahi,b,y1,0,y2,1)+...
    py1x2y2x3(r,s,alo,ahi,b,y1,1,y2,0)+...
    py1x2y2x3(r,s,alo,ahi,b,y1,1,y2,1));

py3x4gy1x2y2x3=@(r,s,alo,ahi,b,y3,x4,y1,x2,y2,x3)TT(r,s,alo,ahi,b,x3,y2,x4,y3); % 4x4 transition probability    

px2x3gy1y2=@(r,s,alo,ahi,b,x2,x3,y1,y2)py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,x3)./py1y2(r,s,alo,ahi,b,y1,y2); % Bayes' rule

py3x4gy1y2=@(r,s,alo,ahi,b,y3,x4,y1,y2)(... % sum over x2,x3=0,1
    py3x4gy1x2y2x3(r,s,alo,ahi,b,1,x4,y1,0,y2,0).*px2x3gy1y2(r,s,alo,ahi,b,0,0,y1,y2)+...
    py3x4gy1x2y2x3(r,s,alo,ahi,b,1,x4,y1,0,y2,1).*px2x3gy1y2(r,s,alo,ahi,b,0,1,y1,y2)+...
    py3x4gy1x2y2x3(r,s,alo,ahi,b,1,x4,y1,1,y2,0).*px2x3gy1y2(r,s,alo,ahi,b,1,0,y1,y2)+...
    py3x4gy1x2y2x3(r,s,alo,ahi,b,1,x4,y1,1,y2,1).*px2x3gy1y2(r,s,alo,ahi,b,1,1,y1,y2)); 

py3gy1y2=@(r,s,alo,ahi,b,y3,y1,y2)(... % sum over x4=0,1
    py3x4gy1y2(r,s,alo,ahi,b,y3,0,y1,y2)+...
    py3x4gy1y2(r,s,alo,ahi,b,y3,1,y1,y2));


HYUB2=@(r,s,alo,ahi,b)(...
    py1y2(r,s,alo,ahi,b,0,0).*phi2(py3gy1y2(r,s,alo,ahi,b,1,0,0))+...
    py1y2(r,s,alo,ahi,b,0,1).*phi2(py3gy1y2(r,s,alo,ahi,b,1,0,1))+...
    py1y2(r,s,alo,ahi,b,1,0).*phi2(py3gy1y2(r,s,alo,ahi,b,1,1,0))+...
    py1y2(r,s,alo,ahi,b,1,1).*phi2(py3gy1y2(r,s,alo,ahi,b,1,1,1)));
MIUB2=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYUB2(r,s,alo,ahi,b);
MIUB2plot=min(MIUB2(rplot,splot,aloplot,ahiplot,bplot),HX(rplot,splot));
figure
[c,h]=contour(rplot,splot,MIUB2plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Upper Bound 2: min(H(X)-H(X,Y)+H(Y3|Y2,Y1),H(X))    ','FontSize',20)
shg
%print -dpdf Cap2StateMarkov2_UB2.pdf

%% Next simplest Lower Bound is based on H[Y3|X2,Y1,Y2]

py1x2y2=@(r,s,alo,ahi,b,y1,x2,y2)(... % sum over x3
    py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,0)+...
    py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,1));

px3gy1x2y2=@(r,s,alo,ahi,b,x3,y2,x2,y1)(...
    py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,x3)./py1x2y2(r,s,alo,ahi,b,y1,x2,y2));

py3x4gx2y1y2=@(r,s,alo,ahi,b,y3,x4,x2,y1,y2)(...% sum over x3=0,1
    py3x4gy1x2y2x3(r,s,alo,ahi,b,y3,x4,y1,x2,y2,0).*px3gy1x2y2(r,s,alo,ahi,b,0,y2,x2,y1)+...
    py3x4gy1x2y2x3(r,s,alo,ahi,b,y3,x4,y1,x2,y2,1).*px3gy1x2y2(r,s,alo,ahi,b,1,y2,x2,y1));

px2y1y2=@(r,s,alo,ahi,b,x2,y1,y2)(...% sum over x3=0,1
    py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,0)+...
    py1x2y2x3(r,s,alo,ahi,b,y1,x2,y2,1));

py3gx2y1y2=@(r,s,alo,ahi,b,y3,x2,y1,y2)(...% sum over x4=0,1
    py3x4gx2y1y2(r,s,alo,ahi,b,y3,0,x2,y1,y2)+...
    py3x4gx2y1y2(r,s,alo,ahi,b,y3,1,x2,y1,y2));

HYLB2=@(r,s,alo,ahi,b)(... % sum over x2,y1,y2=0,1
    px2y1y2(r,s,alo,ahi,b,0,0,0).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,0,0,0))+...
    px2y1y2(r,s,alo,ahi,b,0,0,1).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,0,0,1))+...
    px2y1y2(r,s,alo,ahi,b,0,1,0).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,0,1,0))+...
    px2y1y2(r,s,alo,ahi,b,0,1,1).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,0,1,1))+...
    px2y1y2(r,s,alo,ahi,b,1,0,0).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,1,0,0))+...
    px2y1y2(r,s,alo,ahi,b,1,0,1).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,1,0,1))+...
    px2y1y2(r,s,alo,ahi,b,1,1,0).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,1,1,0))+...
    px2y1y2(r,s,alo,ahi,b,1,1,1).*phi2(py3gx2y1y2(r,s,alo,ahi,b,1,1,1,1)));
MILB2=@(r,s,alo,ahi,b)HX(r,s)-HXY(r,s,alo,ahi,b)+HYLB2(r,s,alo,ahi,b);
MILB2plot=MILB2(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,MILB2plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Lower Bound 2: H(X)-H(X,Y)+H(Y3|X2,Y1,Y2)','FontSize',20)
shg

%% Compare LB2 and UB2

bound_diff_2=@(r,s,alo,ahi,b)min(HX(r,s),MIUB2(r,s,alo,ahi,b))-MILB2(r,s,alo,ahi,b);
bd2plot=bound_diff_2(rplot,splot,aloplot,ahiplot,bplot);
figure
[c,h]=contour(rplot,splot,bd2plot);
clabel(c,h);
axis equal,axis tight,axis([0 1 0 1])
colorbar
set(gca,'FontSize',20)
title('Difference between UB2 and LB2  ','FontSize',20)
shg

%% compare contour lines of LB2 and UB2
figure
[cl,hl]=contour(rplot,splot,MILB2plot,'r--');
clabel(cl,hl);
hold on
[cu,hu]=contour(rplot,splot,MIUB2plot,'k-');
clabel(cu,hu);
% add line along IID axis
line([0 1],[1 0])
axis equal, axis tight, axis([0 1 0 1])
set(gca,'FontSize',20)
title('Lower and Upper Bounds of MI (Depth=2)  ','FontSize',20)
shg
print -dpdf Capacity2StateMarkov2_alo_p1_ahi_p9_b_b5.pdf

%% Find maximum of lower bound (roughly) over r and s grid
LB2max_rough=max(max(MILB2(rplot,splot,aloplot,ahiplot,bplot))); % 0.278897344739580
UB2min_rough=min(min(min(HX(rplot,splot),MIUB2(rplot,splot,aloplot,ahiplot,bplot))));
UB2max_rough=max(max(min(HX(rplot,splot),MIUB2(rplot,splot,aloplot,ahiplot,bplot)))); % 0.278897850764829
% difference is 5.0603e-07
figure
%[c,h]=contour(rplot,splot,MIUB2(rplot,splot,aloplot,ahiplot,bplot),[LB2max_rough,LB2max_rough]);
surf(rplot,splot,MILB2(rplot,splot,aloplot,ahiplot,bplot));
hold on
surf(rplot,splot,MIUB2(rplot,splot,aloplot,ahiplot,bplot));
shg
shading interp
