% Behavior of information rate for 2-state input, 2-state output Markov
% channel as we decrease the time step. The idea is that the
% binding/unbinding rates per time step are proportional to dt.  

%% function definitions
abar=@(alo,ahi,x)(1 - x).*alo + x.*ahi;
H=@(p)p.*log2(1./p)+(1-p).*log2(1./(1 - p));
MIrate=@(alo,ahi,b,x,epsilon)...
    (H(abar(alo.*epsilon,ahi.*epsilon,x))...
    -(x.*H(ahi.*epsilon)+(1-x).*H(alo.*epsilon)))...
    ./...
    (1+(abar(alo.*epsilon,ahi.*epsilon,x)./(b*epsilon)));

MIratescaled=@(alo,ahi,b,x,epsilon)...
    (1./epsilon).* ...
    ((H(abar(alo.*epsilon,ahi.*epsilon,x))...
    - (x.*H(ahi.*epsilon)+(1-x).*H(alo.*epsilon)))...
    ./...
    (1+(abar(alo.*epsilon,ahi.*epsilon,x)./(b*epsilon))));

xplot=linspace(1e-6,1-1e-6,201);

%% first plot (unscaled)
figure
for epsilon=10.^(0:-1:-4)
    plot(xplot,log(MIrate(.1,.9,.5,xplot,epsilon)),'LineWidth',3)
    hold on
end
grid on
axis([0 1 -20 0])
set(gca,'FontSize',20)
xlabel('p_H ','FontSize',20)
ylabel('log ( MI rate ) ','FontSize',20)
shg

print -dpdf Capacity2StateMarkov-scaling-time-p1p9p5.pdf
    
%% second plot (scaled)

figure
for epsilon=10.^(0:-1:-4)
    plot(xplot,MIratescaled(.1,.9,.5,xplot,epsilon),'LineWidth',3)
    hold on
end
grid on
axis([0 1 0 .3])
set(gca,'FontSize',20)
xlabel('p_H ','FontSize',20)
ylabel('MI rate / \epsilon ','FontSize',20)
shg
print -dpdf Capacity2StateMarkov-scaling-time-p1p9p5-scaled.pdf