function [value,isterminal,direction] = impact(t,y)
% Locate the time when tip touches sample in a decreasing direction and stop integration.
global b o

l = b + y(1) - o;           % Distance between sample and tip [m]

value = l;                  % detect height = 0
isterminal = 1;             % stop the integration
direction = -1;             % negative direction
end