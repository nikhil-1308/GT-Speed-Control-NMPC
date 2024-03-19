function [c, ceq] = ConstraintFCN(u,uold,x,Ts,Np,LBo,UBo,LBdu,UBdu,p)

%% Constraint function of nonlinear MPC for F8 system

%

% Inputs:

% u: optimization variable, from time k to time k+Np-1

% x: current state at time k

% Ts: controller sample time

% Np: prediction horizon length

% uold: latest applied control input

% LBo: Lower bound of output x

% UBo: Upper bound of output x

% LBdu: Lower bound for input difference uk - uk-1

% UBdu: Upper bound for input difference uk - uk-1

% p: Parameters for model

%

% Output:

% c: inequality constraints applied across prediction horizon

% ceq: equality constraints (empty)

%

%% Nonlinear MPC design parameters

% range of angle of attack

zMin = LBo;

zMax = UBo;

%% Inequality constraints calculation

c = zeros(2*Np,1);

cy = zeros(2*Np,1);

% Apply 2N constraints across prediction horizon, from time k+1 to k+Np

xk = x;

uk = u(1);

duk = u(1)-uold;

for ct=1:Np

% Obtain new state at next prediction step

% xk1 = GTStateFcnDT(xk,uk);
xk1 = state_sim(xk,uk);

% yk1 = GTOutputFcn(xk1,uk);

% -z + zMin < 0 % lower bound

% cy(2*ct-1) = -yk1(1)+zMin;

% cy(2*ct-1) = -xk1(1)+zMin;
cy(2*ct-1) = -xk1(2)+zMin;

c(2*ct-1) = -duk+LBdu;

% z - zMax < 0 % upper bound

% cy(2*ct) = yk1(1)-zMax;

% cy(2*ct) = xk1(1)-zMax;
cy(2*ct) = xk1(2)-zMax;

c(2*ct) = duk-UBdu;

% Update plant state and input for next step

xk = xk1;

if ct<Np

uk = u(ct+1);

duk = u(ct+1)-u(ct);

end

end

c = [c;cy];

%% No equality constraints

ceq = [];