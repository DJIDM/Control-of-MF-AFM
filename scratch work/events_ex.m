function [value,isterminal,direction] = events_ex(t,y,l)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y(1)-l;     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
