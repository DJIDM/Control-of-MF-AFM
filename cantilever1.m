function dydt = cantilever1(t,y,A)
% For unimodal AFM operation

global D Q1 Ar w1 wd

% Interpolate u to obtain the value of the time-dependent terms at the specified time.
u = D*sin(wd*t);

% Interaction Force
Fts = DMT(y(1));

% State-space equations for impact oscillator AFM cantilever model
dydt = [y(2);                                          ... % Tip velocity
        -((w1)^2)*y(1) - (w1/Q1)*y(2) + u + Fts;         ... % Tip acceleration
        Ar - A;                                        ... % Amplitude error 
       ];
end