function [dx, y] = cantilever_NLGB(t, x, u, b, o, r, w_n, Q, varargin)
% AFM cantilever state-space model.
% t = time [s], x = state vector, u = input signal, 
% b = cantilever base height [m]
% o = sigma = height of sample surface [m]
% r = restitution coefficient []
% w_n = natural frequency [rad/s]
% Q = quality factor []



% Output equation.
y = x(1);                                               % Vertical position [m] of the cantilever tip w/ respect to 'b'.


% Interaction force 
H = 1.4e-19;                                            % Hamaker constant [J]
r_t = 2e-9;                                             % tip radius [m]
l_m = 0.42e-9;                                          % intermolecular distance [m]
E_t = 1.65e11;                                          % elastic modulus (tip) [Pa]
E_s = 1.65e11;                                          % elastic modulus (sample) [Pa]
V_t = 0.27;                                             % Poisson ratio (tip) []
V_s = 0.27;                                             % Poisson ratio (sample) []


l = b + x(1) - o;                                       % Distance between sample and tip [m]

% DMT model for interaction force F
if l > l_m
    F = -(H*r_t)/(6*l^2);
else
    F = -(H*r_t)/(6*l_m^2) + (4/3)*sqrt(r_t*(l_m - l)^3)/((1-V_t^2)/E_t +(1-V_s^2)/E_s);
end


% Reset law
if l == 0
    x(2) = -r*x(2);
end


% State equations.
dx = [x(2);                                         ... % Tip velocity
      -((w_n)^2)*x(1) - (w_n/Q)*x(2) + u(1) + F     ... % Tip acceleration
     ];
 