 

% Automatic flight control of F8 aircraft using model predictive control

% and true system dynamics as model

clear all, close all, clc

warning('off')

global N TauL Ain alpha Cf phi mdotF

load 'XREF'

% initial conditions

Ain = 14.67;

Cf = 350;

alpha = 10;

N = 3200;

TauL = 300; % Load torque

phi = 0.55;

%% FSR Calc

Cv = 0.5;

Sg = 0.6;

dp = 1;

Tref = 85.3549;

FF_FSR = 70;

FF_FSRT = 70;

%% MPC

dt = 0.08;

Ts = dt;

Nvar = 2;

% Parameters MPC

options = optimoptions('fmincon','Algorithm','sqp','Display','none', 'MaxIterations',100);

Duration = 4; % Run for 'Duration' time units

Ton = 0; % Time units when control is turned on

Np = 10;

Nu = Np;

% Q = 0.1;

Q = diag([0.01 0.01]);

R = 0.1;

Ru = 0.1;

LBo = x1(2) - 2e4;

UBo = x1(2);

LBdu = 0.2;

UBdu = 0.75;

LB = LBdu;

UB = UBdu;

x0n=[0 16800]'; % Initial condition

Results = struct('r', [], 'x',[], 'u', [], 't', [], 'xref', [], 'J', []);

%% Run MPC

% Prepare variables // Initialization

Nt = (Duration/Ts)+1;

tspan = 0:Ts:Duration;

uopt0 = phi;

Mdot0 = 0.0041;

xhat = x0n;

uopt = uopt0.*ones(Nu,1);MdotfHistory = Mdot0.*ones(Nu,1);

xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;

uHistory = zeros(1,Nt); uHistory(1) = uopt(1);

tHistory = zeros(1,Nt); tHistory(1) = 0;

rHistory = zeros(2,Nt); yHistory = zeros(1,Nt); %rHistory = zeros(1,Nt);

% Start simulation

fprintf('Simulation started. It might take a while...\n')

tic

for ct = 1:Duration/Ts

% Set references over optimization horizon

tref = (ct:ct+Np-1).*Ts;

xref = x1;

if tspan(ct) >= 3+dt

    % NMPC with full-state feedback

    COSTFUN = @(u) ObjectiveFCN(u,xhat,Np,xref,uHistory(:,ct),[],Q,R,Ru);

    CONSFUN = @(u) ConstraintFCN(u,uHistory(:,ct),xhat,Ts,Np,LBo,UBo,LBdu,UBdu,[]);

    uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,options);

end

if tspan(ct) >= 3+2*dt

    TauL = 310;

elseif tspan(ct) >= 3

    TauL = 190;

end

% Integrate system

% xhat = GTStateFcnDT(xhat,uopt(1)); % Increase time resolution for simulation & keep control constant

xhat = state_sim(xhat,uopt(1)); % Increase time resolution for simulation & keep control constant

yhat = GTOutputFcn(xhat,uopt(1));

 

if ct == 1

    Anm = (((mdotF/Cv)^2)*Sg)/dp;

end

 

[QF,Fsr,Fsrt] = SimFSR(Tref,FF_FSR,FF_FSRT,Cv,Sg,dp,Anm);

 

FF_FSR = Fsr;

FF_FSRT = Fsrt;

Mf(:,ct+1) = mdotF;

MdotfHistory(:,ct+1) = QF;

xHistory(:,ct+1) = xhat;

uHistory(:,ct+1) = uopt(1);

tHistory(:,ct+1) = ct*Ts;

rHistory(:,ct+1) = xref;

yHistory(:,ct+1) = yhat;

end

tElapsed = toc;

fprintf('Simulation finished!\n')

%% Collect results

Results.eval_time = tElapsed;

Results.xref = rHistory;

Results.x = xHistory;

Results.u = uHistory;

Results.t = tHistory;

% Results.J = evalObjectiveFCN(Results.u,Results.x,Results.xref,diag(Q),R,Ru);

%% Show results

figure,plot(xHistory')%, hold on, plot(r,'-k'), legend('aoa','pa','pr','r','y')

figure,plot(uHistory)

% figure,plot(Results.J)

figure,plot(yHistory)