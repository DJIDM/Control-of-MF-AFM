function F = DMT(y)
% Derjaguin-Muller-Toporov model for interaction forces in AFM.
% y = tip height [m]
% b = cantilever base height [m]
% o = sigma = height of sample surface [m]

global b o

% Parameters 
H = 1.4e-19;                            % Hamaker constant [J]
r_t = 2e-9;                             % tip radius [m]
l_m = 0.42e-9;                          % intermolecular distance [m]
E_t = 1.70e11;                          % elastic modulus (tip) [Pa]
E_s = 1.65e11;                          % elastic modulus (sample) [Pa]
V_t = 0.27;                             % Poisson ratio (tip) []
V_s = 0.27;                             % Poisson ratio (sample) []


l = b + y - o;                          % Distance between sample and tip [m]

% DMT model for interaction force F
if l > l_m
    F = -(H*r_t)/(6*l.^2);
else
    F = -(H*r_t)/(6*l_m^2) + (4/3)*sqrt(r_t*(l_m - l).^3)/((1-V_t^2)/E_t +(1-V_s^2)/E_s);
end