% ------------------------------ INITIALIZATION --------------------------------------
clf; clear; % close all
addpath 'C:\Users\d3260\OneDrive - University of Bristol\Research Skills\TB2\Research Project\Corragio'

% -------------------------------- USER INPUT ----------------------------------------
fprintf('Sample file to load '); Nf = input('(1 = DNA, 2 = TiS2, 3 = UO): ');
switch Nf
    case 1
        load('DNA.mat')
    case 2
        load('titaniumDisulfide.mat')
    case 3
        load('uraniumOxide.mat')
    otherwise
        error('File DNE')
end

Nr = size(trueSample,1);                               % # of rows of selected sample
fprintf('Line no. to be sampled (1 - %u)',Nr); Nl = round(input(': '));
if Nl < 1 || Nl > Nr
    error('Line DNE')
end

x0 = input('Initial state of cantilever ([pos (m); vel (m/s)]): ');
dt = input('Sampling interval for input signal (sec): ');
tT = input('Duration of input signal (sec): '); 
A_f = input('Free oscillation amplitude of cantilever (m): ');
w_d = input('Driving frequency of dither piezo (rad/s): ');


% --------------------------- CANTILEVER PARAMETERS ----------------------------------
Q = 100;                                                % quality factor []
k = 42;                                                 % stiffness coefficient [N/m]
m = 1.3098e-11;                                         % mass [kg]
r = 0.9;                                                % restitution coefficient

w_n = sqrt(k/m);                                        % natural frequency [rad/s]
b = pid(0,10000,0);                                     % 'b' is estimated as a PID controller (K_P = 0, K_I = 10^4, K_D = 0)

% ---------------------------- DITHER PIEZO SIGNAL -----------------------------------
o = trueSample(Nl,:);                                   % sigma from cantilever.m = user chosen line no. of loaded sample
A_r = 0.9*A_f;                                          % amplitude reference value
D = A_f*abs(w_n^2 - w_d^2 + (w_n/Q)*(1i*w_d));          % driving amplitude of dither piezo input signal
T_d = (2*pi)/w_d;                                       % period of input signal
[u,t] = gensig('sin',T_d,tT,dt);                        % generate input signal 
u = iddata([],D*u,dt,'Name','Dither Piezo Action');     % set input signal amplitude to 'D'
plot(u)                                                 % visualize input signal


% --------------------------------- AFM Model ----------------------------------------
% Create model
AFM_sys = idnlgrey('cantilever_NLGB',[1 1 2],[(A_r - A_f); o(1); r; w_n; Q], x0, 0, ...
'Name','AFM', ...
'InputName','Dither Piezo Action','InputUnit', 'm/(s^2)', ...
'OutputName', 'Cantilever Tip Vertical Position','OutputUnit', 'm', ...
'TimeUnit', 's');
AFM_sys = setinit(AFM_sys, 'Name', {'Cantilever Tip Vertical Position', ...
                                    'Cantilever Tip Vertical Velocity'});
AFM_sys = setinit(AFM_sys, 'Unit', {'m', 'm/s'});
AFM_sys = setpar(AFM_sys, 'Name', {'b: Cantilever Base Height', ...
                                   'o: Sample Surface Height', ...
                                   'r: Restitution Coefficient', ...
                                   'w_n: Natural Frequency', ...
                                   'Q: Quality Factor'});
                               
% Set r, w_n, and Q as fixed parameters
for i = 3:5   
    AFM_sys.Parameters(i).Fixed = true;
end

present(AFM_sys);                                         % View model properties