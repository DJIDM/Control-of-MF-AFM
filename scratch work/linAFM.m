% ------------------------------ INITIALIZATION ---------------------------------
clf; clear; % close all

% -------------------------------- PARAMETERS -----------------------------------
% Cantilever
w_n = 2.85e5*2*pi;            % natural frequency [rad/s]
Q = 100;                      % quality factor []
r = 0.9;                      % impact coefficient of restitution []
k = 42;                       % stiffness coefficient [N/m]
m = k/(w_n)^2;                % mass [kg]
c = m*w_n/Q; 


% Forces 
H = 1.4e-19;                  % Hamaker constant [J]
r_t = 2e-9;                   % tip radius [m]
l_m = 0.42e-9;                % intermolecular distance [m]
E_t = 1.65e11;                % elastic modulus (tip) [Pa]
E_s = 1.65e11;                % elastic modulus (sample) [Pa]
V_t = 0.27;                   % Poisson ratio (tip) []
V_s = 0.27;                   % Poisson ratio (sample) []



% ----------------------------------- MODEL --------------------------------------
% DMT
l = 0:1e-12:3e-9; F = zeros(size(l));
for i = 1:length(l)
    if l(i) > l_m
        F(i) = -(H*r_t)/(6*l(i)^2);
    else
        F(i) = -(H*r_t)/(6*l_m^2) + (4/3)*sqrt(r_t*(l_m - l(i))^3)/((1-V_t^2)/E_t +(1-V_s^2)/E_s);
    end
end

% Plotting the interaction force 'F' vs tip-sample seperation 'l'
plot(l*1e9, F*1e9);
title( {'Interaction forces (DMT model)'}, 'Interpreter','latex', 'FontSize', 10 );
xlabel( {'Tip-sample distance $l$ [nm]'}, 'Interpreter','latex', 'FontSize', 10 );
ylabel( {'Force $F$ [nN]'}, 'Interpreter','latex', 'FontSize', 10 );
grid on;


% State-space representation
cantileverA = [0,1;-(w_n)^2,-(w_n)/Q];
cantileverB = [0;1];
cantileverC = [1,0];
cantileverD = [0];
cantilever = ss(cantileverA,cantileverB,cantileverC,cantileverD);

% Simulation
figure
lsim(cantilever)                      % dither piezo response
