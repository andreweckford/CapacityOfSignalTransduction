% This MATLAB script generates a Monte Carlo mutual information estimate
% for every point on a 2-d grid representing values of r and s.
% 
% The default input process parameters will produce an 79x79 matrix of 
% mutual information estimates, where the (i,j)th entry corresponds
% to input process transition probabilities r = i/80 and s = j/80.
% The default quality parameters give a HIGH QUALITY estimate of 
% mutual information, though with VERY LONG running time (several days).
%
% At the end of this routine, the resulting matrix is stored in the
% variable z -- see below
%
% Written by Andrew W. Eckford, aeckford@yorku.ca

clear

%%%%% Parameters

%% Model parameters
%
% binding probability at low concentration
alphaL = 0.1; 

% binding probability at high concentration
alphaH = 0.5; 

% unbinding probability
beta = 0.2;

%% Input process parameters
%
% r and s are the transition probabilities for the input process
% num contains the number of points per dimension in the 2-dimensional
% grid representing r and s
num = 79; 

% step is the size of the quantization interval in r and s
% r and s run from step to num*step
% thus, step should normally be 1/(num+1)
% moreover, step should NOT be greater than 1/(num+1)
step = 1/80; 

%% Quality parameters -- increase for higher quality (but longer run time)
%
% length of the vector of channel uses used in each iteration of Monte Carlo
mcLength = 10000;

% number of iterations of Monte Carlo per point on the r-s grid
mcIter = 1000;

%%%%% End of parameters

PB = zeros(2,2,2);
PB(:,:,1) = [1-alphaL alphaL; beta 1-beta];
PB(:,:,2) = [1-alphaH alphaH; beta 1-beta];
z = zeros(num);
rr = step:step:(num*step);
ss = step:step:(num*step);
for c=1:num
    for d = 1:num
        r = rr(c);
        s = ss(d);
        P = [1-r r; s 1-s];
        z(c,d) = mutualInfoMC(P,PB,[1 0],mcLength,mcIter);
    end
end

% the result is stored in z
