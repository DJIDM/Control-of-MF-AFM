function [b,u] = hPID(A,Af,Ar,tA,t,wn,wd,v,u,Q)
% Hybrid PID controller
% b = new base height; u = new dither piezo input signal
% A = tip oscillation amplitude [m] ~ (base height (b) - sample height (o))
% Af = free oscillation amplitude [m]
% Ar = reference amplitude [m]
% tA = time constant [s]
% t = time vector [s]
% wn = cantilever natural frequency [rad/s]
% wd = dither piezo driving frequency [rad/s]
% v = tip velocity (y2) [m/s]
% u = input signal
% Q = quality factor

last = @(V) V(end);                     % function to retrieve last entry in an vector

%% PARAMETERS -------------------------------------------------------------
KP = 0;                                 % proportional gain
KI = 10000;                             % integral gain
KD = 0;                                 % derivative gain
Kt = 2.3;                               % time constant gain
Q_ = 30;                                % Q-control effective quality factor
dQPL = 25;                              % probe loss Q
dQR = 25;                               % recoil Q
At = 0.95*Af;                           % threshold amplitude (high)
At_ = 0.94*Af;                          % threshold amplitude (low)
AtRL = 0.5*Af;                          % threshold amplitude - recoil
a = -400*Af;                            % signal noise threshold
q = 1;                                  % start in Regular mode

%% PROBELOSS & RECOIL Q ---------------------------------------------------
QPL = Q_ - dQPL*min([abs((Ar - A)/(Ar - Af)),1]);
QR = Q_ - dQR*min([abs((Ar - A)/(Ar - Af)),1]);

%% PID CONTROL ------------------------------------------------------------
PID = @(z) KP*z(end) + KI*trapz(1/wd,z) + KD*last(gradient(z));

%% SCHEME -----------------------------------------------------------------
switch q
    case 1                                                   % Regular mode
        Ks = 1;                                              % dynamic PID gain
        D_ = Af*abs(wn^2 + (wn/Q_)*(1i*wd) + (1i*wd)^2);     % alternate input amplitude
        KQ = wn*(1/Q_ - 1/Q);                                % Q-control gain
    case 2                                                   % ProbeLoss mode
        Ks = 15;
        D_ = Af*abs(wn^2 + (wn/Q_)*(1i*wd) + (1i*wd)^2);
        KQ = wn*(1/Q_ - 1/Q);
    case 3                                                   % Recovery mode
        Ks = 1;
        D_ = Af*abs(wn^2 + (wn/QPL)*(1i*wd) + (1i*wd)^2);
        KQ = wn*(1/QPL - 1/Q);
    case 4                                                   % Recoil mode
        Ks = 1;
        D_ = Af*abs(wn^2 + (wn/QR)*(1i*wd) + (1i*wd)^2);
        KQ = wn*(1/QR - 1/Q);
end

%% TRANSITION GUARDS ------------------------------------------------------
dA = last(gradient(A));                                    % A dot

if q == 1 
    if A >= At                                             % g1,2
        q = 2;
    else
    if A <= AtRL                                           % g1,4
        q = 4;
        p = 0;
    end
    end
end

if q == 2
    if A <= At_                                            % g2,1
        q = 1;
    else
    if dA < a                                              % g2,3
        q = 3;
        p = 0;
    end
    end
end

if q == 3
    if (dA > 0) && ~p
        q = 3;
        p = 1;
    else
    if ((dA < 0) && p) || ((t(end) - t(1)) >= Kt*tA)       % g3,1
        q = 1;
    end
    end
end

if q == 4
    if (dA > 0) && ~p                                      % g4,4
        q = 4;
        p = 1;
    else
    if ((dA < a) && p) || ((t(end) - t(1)) >= Kt*tA)       % g4,1
        q = 1;
    end
    end
end

%% OUTPUT -----------------------------------------------------------------
% Dynamic PID
if A <= At
    b = PID(Ar - A);
else
    b = PID((Ar - At) + Ks*(At - A));
end

% Q-control
u = D_*u - KQ*v;

end