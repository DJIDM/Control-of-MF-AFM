%% ------------------------------ INITIALIZATION --------------------------------------
clear; % close all

global m b o k1 k2 w1 w2 wd D F1 F2 Q1 Q2 Ar

last = @(V) V(end);                                            % function to retrieve last entry in an vector

% cantilever
m = 1.3098e-11;                                                % mass [kg]
c = 2.3455e-7;                                                 % damping coefficient [kg/s]
r = 0.9;                                                       % restitution coefficient

Nm = input('Operation Mode (1 = monomodal, 2 = bimodal): ');
switch Nm
    case 1
        %% -------------------------------- USER INPUT 1 ----------------------------------------
fprintf('Sample file to load '); Nf = input('(0 = Flat, 1 = Ramp, 2 = Wave, 3 = DNA, 4 = TiS2, 5 = UO): ');
switch Nf
    case 0
        load('flatSurface.mat')
        Af = 50e-9;
    case 1
        load('rampSurface.mat')
        Af = 50e-9;
    case 2
        load('sqWave.mat')
        Af = 50e-9;
    case 3
        load('DNA.mat')
        Af = 2e-9;
    case 4
        load('titaniumDisulfide.mat')
        Af = 2e-9;
    case 5
        load('uraniumOxide.mat')
        Af = 200e-9;
    otherwise
        error('File DNE')
end

Nc = size(trueSample,2);                                       % # of columns of selected sample
fprintf('Line no. to be sampled (1 - %u)',Nc); Nl = round(input(': '));
if Nl < 1 || Nl > Nc
    error('Line DNE')
end

tT = input('Duration of input signal (sec): ');
wd = input('Driving frequency of dither piezo (rad/s): ');

%% -------------------------------- PARAMETERS 1 ---------------------------------------
% cantilever
k1 = 42;                                                       % stiffness coefficient [N/m]
w1 = sqrt(k1/m);                                               % natural frequency [rad/s]
Q1 = m*w1/c;                                                   % quality factor []

%% -------------------------------- DITHER PIEZO 1 --------------------------------------
Ar = 0.9*Af;                                                   % amplitude reference value
D = Af*abs(w1^2 - wd^2 + (w1/Q1)*(1i*wd));                     % driving amplitude of dither piezo input signal
% tA = (2*Q1)/w1;                                              % dither piezo time constant

%% -------------------------------- Z-AXIS PIEZO 1 --------------------------------------
b = [Ar];                                                      % initial base height
figure; mesh(trueSample)
s = trueSample(:,Nl);                                          % line of sample to be scanned chosen by user
figure; h = linspace(0,length(s),length(s)); plot(h,s); xlim([0 length(s)]);
xlabel('distance [m]'); ylabel('height [m]');
title('Sample Topography');
KP = 0.5; KI = 8.5; KD = 0;                                    % PID controller gains                      

%% ------------------------------ INITIAL CONDITIONS 1 ---------------------------------------
A = b - s(1);
t0 = 0;
y0 = [Ar;0;Ar - A];
refine = 4;
options = odeset('Events',@impact,'OutputFcn',@odeplot,'OutputSel',1,'Refine',refine);

diagram = figure;
axis([0,tT,-inf,inf]);
box on; hold on;

tt = t0; yt = y0.'; tet = []; yet = []; iet = [];       % initialize ODE parameters
eAt = []; bt = [b]; st = [s(1)];     

%% ---------------------------------- SIMULATION 1 ---------------------------------------
for i = 1:length(s)
    o = s(i); 
   % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode15s(@(t,y) cantilever1(t,y,A),[t0 tT],y0,options);
   if ~ishold
      hold on
   end
   % Accumulate output.  This could be passed out as output arguments.
   c = length(t);
   tt = [tt; t(2:c)];
   yt = [yt; y(2:c,:)];
   tet = [tet; te];          
   yet = [yet; ye];
   iet = [iet; ie];
   eAt = [eAt; Ar - A];
   
   % PID controller for base height adjustment
   b = KP*(Ar - A) + KI*y(c,3) + KD*last(gradient(eAt));
   bt = [bt;b];
   A = b - o;
   
   if diagram.UserData.stop
      break;
   end
   
   y0(1) = y(c,1);   y0(2) = -r*y(c,2);                     % Reset Law
   
   
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'Events',@impact,'InitialStep',t(c)-t(c-refine),'MaxStep',t(c)-t(1));
   
   t0 = t(c);
   if t0 == tT
       break
   end
end

plot(tet,yet(:,1),'ro')
xlabel('time');
ylabel('height');
title('Cantilever Tip Trajectory');
hold off
odeplot([],[],'done');
figure; plot(linspace(0,tT,length(eAt)),eAt); title('Error Signal vs. Time'); xlabel('Time (s)'); ylabel('Amplitude Error (m)');xlim([0 tT]);

    case 2
%% -------------------------------- USER INPUT 2 ----------------------------------------
fprintf('Sample file to load '); Nf = input('(0 = Flat, 1 = Ramp, 2 = Wave, 3 = DNA, 4 = TiS2, 5 = UO): ');
switch Nf
    case 0
        load('flatSurface.mat')
        A01 = 10e-9;
    case 1
        load('rampSurface.mat')
        A01 = 6e-9;
    case 2
        load('sqWave.mat')
        A01 = 50e-9;
    case 3
        load('DNA.mat')
        A01 = 2e-9;
    case 4
        load('titaniumDisulfide.mat')
        A01 = 2e-9;
    case 5
        load('uraniumOxide.mat')
        A01 = 200e-9;
    otherwise
        error('File DNE')
end

Nc = size(trueSample,2);                                        % # of columns of selected sample
fprintf('Line no. to be sampled (1 - %u)',Nc); Nl = round(input(': '));
if Nl < 1 || Nl > Nc
    error('Line DNE')
end

tT = input('Duration of input signal (sec): ');
p = input('Free amplitude of 2nd mode (% of 1st mode free amplitude): ');

%% -------------------------------- PARAMETERS 2 ---------------------------------------
% cantilever
A02 = (p/100)*A01;                                              % free amplitude for 2nd mode [m]
f1 = 48.913; f2 = 306.194;                                      % resonant frequencies [kHz]
w1 = 2*pi*f1; w2 = 2*pi*f2;                                     % angular frequencies [rad/s]
k1 = 0.9; k2 = 35.2;                                            % stiffness coefficients [N/m]
Q1 = 255; Q2 = 1000;                                            % quality factor []

%% -------------------------------- DITHER PIEZO 2 --------------------------------------
Ar = 0.9*A01;                                                   % amplitude reference value
F1 = k1*A01/Q1; F2 = k2*A02/Q2;                                 % external excitation forces

%% -------------------------------- Z-AXIS PIEZO 2 --------------------------------------
b = [Ar];                                                       % initial base height
figure; mesh(trueSample)
s = trueSample(:,Nl);                                           % line of sample to be scanned chosen by user
figure; h = linspace(0,s(end),length(s)); plot(h,s)
xlabel('distance [m]'); ylabel('height [m]');
title('Sample Topography');
KP = 0.5; KI = 8.5; KD = 0;                                     % PID controller gains                      

%% ------------------------------ INITIAL CONDITIONS 2 ---------------------------------------
A = b - s(1);
t0 = 0;
y0 = [Ar;0;(0/100)*Ar;0;Ar - A];
refine = 4;
options = odeset('Events',@impact,'OutputFcn',@odeplot,'OutputSel',1,'Refine',refine);

diagram = figure;
axis([0,tT,-inf,inf]);
box on; hold on;

tt = t0; yt = y0.'; tet = []; yet = []; iet = [];               % initialize ODE parameters
eAt = []; bt = [b]; st = [s(1)];     

%% ---------------------------------- SIMULATION 2 ---------------------------------------
for i = 1:length(s)
    o = s(i); 
   % Solve until the first terminal event.
   [t,y,te,ye,ie] = ode15s(@(t,y) cantilever2(t,y,A),[t0 tT],y0,options);
   if ~ishold
      hold on
   end
   % Accumulate output.  This could be passed out as output arguments.
   c = length(t);
   tt = [tt; t(2:c)];
   yt = [yt; y(2:c,:)];
   tet = [tet; te];          
   yet = [yet; ye];
   iet = [iet; ie];
   eAt = [eAt; Ar - A];
   
   % PID controller for base height adjustment
   b = KP*(Ar - A) + KI*y(c,3) + KD*last(gradient(eAt));
   bt = [bt;b];
   A = b - o;
   
   if diagram.UserData.stop
      break;
   end
   
   y0(1) = y(c,1);   y0(2) = -r*y(c,2);                         % Reset Law
   
   
   % A good guess of a valid first timestep is the length of the last valid
   % timestep, so use it for faster computation.  'refine' is 4 by default.
   options = odeset(options,'Events',@impact,'InitialStep',t(c)-t(c-refine),'MaxStep',t(c)-t(1));
   
   t0 = t(c);
   if t0 == tT
       break
   end
end

plot(tet,yet(:,1),'ro')
xlabel('time');
ylabel('height');
title('Cantilever Tip Trajectory');
hold off
odeplot([],[],'done');
figure; plot(tt,yt(:,3),'-o'); xlim([0 tT]);
figure; plot(linspace(0,tT,length(eAt)),eAt); title('Error Signal vs. Time'); xlabel('Time (s)'); ylabel('Amplitude Error (m)');xlim([0 tT]);

otherwise
        error('Mode DNE')
end
