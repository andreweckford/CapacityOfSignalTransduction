function y = bindingLogP(b,P,PB,py,px)

sx = size(P);
numStates = sx(1);

if (nargin == 4)
    foo = P^1000;
    px = foo(1,:);
end

% initial sum-product message
% at the right-going output of the first state variable node

msg = zeros(numStates,2);

for c = 1:numStates
    if (b(1) == 1)
        msg(c,:) = [px(c)*py(1) 0];
    else
        msg(c,:) = [0 px(c)*py(2)];
    end
end

% norm is used to store the sum of the logs of
% message normalizations (to avoid underflow)
norm = 0;

for t = 2:length(b)
    
    % message through factor node to next state node
    nextMsg = zeros(numStates,2);
    
    for prevState = 1:numStates
        for prevBind = 1:2
            for nextState = 1:numStates
                for nextBind = 1:2
                    nextMsg(nextState,nextBind) = ...
                        nextMsg(nextState,nextBind) ...
                        + msg(prevState,prevBind) ...
                        * P(prevState,nextState) ...
                        * PB(prevBind,nextBind,prevState);
                end
            end
        end
    end
    
    msg = nextMsg;
    
    % message through next state node to next factor node
    % we can observe b, but not the state
    
    for c = 1:numStates
        if (b(t) == 1)
            msg(c,:) = [msg(c,1) 0];
        else
            msg(c,:) = [0 msg(c,2)];
        end
    end
    
    % numerical normalization step: 
    % normalize the messages so they sum to 1
    % and store the sum of the log of the normalization constants
    
    norm = norm + log(sum(sum(msg)));
    msg = msg / sum(sum(msg));
    
end
    
%y = log(sum(sum(msg))) + norm;
% actually sum(sum(msg)) = 1 thanks to normalization

 y = norm;