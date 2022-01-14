function dydt = cantilever2(t,y,A)
% For bimodal AFM operation

global m k1 k2 w1 w2 F1 F2 Q1 Q2 Ar

% Interpolate u to obtain the value of the time-dependent terms at the specified time.
u1 = F1*cos(w1*t); u2 = F2*cos(w2*t); 

% Interaction Force
Fts = DMT(y(1)+y(3));

% State-space equations for impact oscillator AFM cantilever model
dydt = [y(2);                                                  ... % 1st mode velocity
        -(k1/m)*y(1) - (w1/Q1)*y(2) + u1/m + u2/m + Fts;       ... % 1st mode acceleration
        y(4);                                                  ... % 2nd mode velocity
        -(k2/m)*y(3) - (w2/Q2)*y(4) + u1/m + u2/m + Fts;       ... % 2nd mode acceleration
        Ar - A;                                                ... % 1st mode amplitude error
       ];
end
