tmax=1e6; % number of iterations for a movie
t=1;
longrun=1; % set to 1 for long run with no visualization

%% (r,s) grid 
nr=11;
ns=10;
[r,s]=meshgrid(linspace(1/(2*nr),1-1/(2*nr),nr),linspace(1/(2*ns),1-1/(2*ns),ns));

%%  initial conditions
x=nan(ns,nr,tmax);
x(:,:,1)=rand(ns,nr)>s./(r+s);
% visualize them
pcolor(r,s,double(x(:,:,1))),shading flat, colormap gray,axis equal,axis tight,shg

if ~longrun
%% generate Markov chain ensemble "X" alone

while t<tmax
    xold=x(:,:,t);
    t=t+1;
    rtmp=rand(ns,nr);
    x(:,:,t)=xold+(1-xold).*(rtmp<r)-xold.*(rtmp<s);
end

%% visualize as a movie

for t=1:tmax
    pcolor(r,s,double(x(:,:,t))),shading flat, colormap gray,axis equal,axis tight
    title(num2str(t),'FontSize',20)
    MovieX(t) = getframe;
end
movie(MovieX)

%% clean up from first movie
clear M x

end % if ~longrun

%% constants for (X,Y) process
alo=.1;
ahi=.9;
b=.5;
abar=(alo*s+ahi*r)./(r+s); % average value of alpha over X at steady state
dalpha=ahi-alo; % difference in alpha

%%  initial conditions, again
t=1;
x=nan(ns,nr,tmax);
y=nan(ns,nr,tmax);
x(:,:,1)=rand(ns,nr)>s./(r+s);
y(:,:,1)=rand(ns,nr)>b./(b+abar);

%% visualize the initial conditions

    subplot(1,2,1)
    pcolor(r,s,double(x(:,:,t))),shading flat, colormap gray,axis equal,axis tight
    title(num2str(t),'FontSize',20)
    xlabel('X','FontSize',20)
    subplot(1,2,2)
    pcolor(r,s,double(y(:,:,t))),shading flat, colormap gray,axis equal,axis tight
    %title(num2str(t),'FontSize',20)
    xlabel('Y','FontSize',20)


%% generate Y along with X
while t<tmax
        
    if ~rem(t,1000)
        disp(num2str(t))
        save Capacity2StateMarkovMonteCarloData.mat
    end

    xold=x(:,:,t);
    yold=y(:,:,t);
    rtmpx=rand(ns,nr); % for X transitions
    rtmpy=rand(ns,nr); % for Y transitions
    x(:,:,t+1)=xold+(1-xold).*(rtmpx<r)-xold.*(rtmpx<s);
    alpha=alo+dalpha*x(:,:,t+1);
    y(:,:,t+1)=yold+(1-yold).*(rtmpy<alpha)-yold.*(rtmpy<b);    
    t=t+1;
end
save Capacity2StateMarkovMonteCarloData.mat
%% visualize as a movie

if ~longrun
for t=1:tmax
    subplot(1,2,1)
    pcolor(r,s,double(x(:,:,t))),shading flat, colormap gray,axis equal,axis tight
    title(num2str(t),'FontSize',20)
    xlabel('X','FontSize',20)
    subplot(1,2,2)
    pcolor(r,s,double(y(:,:,t))),shading flat, colormap gray,axis equal,axis tight
    %title(num2str(t),'FontSize',20)
    xlabel('Y','FontSize',20)
    MovieXY(t) = getframe;
end
movie(MovieXY)
end % if ~longrun
