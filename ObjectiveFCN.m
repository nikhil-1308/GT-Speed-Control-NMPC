function J = ObjectiveFCN(u,x,Np,xref,u0,p,Q,R,Ru)

%% Cost function of nonlinear MPC for Lotka-Volterra system

%

% Inputs:

% u: optimization variable, from time k to time k+Np-1

% x: current state at time k

% Ts: controller sample time

% Np: prediction horizon

% Nu: control horizon [not implemented]

% xref: state references, varying from time k+1 to k+Np

% u0: previous controller output at time k-1

% p: Parameters for model

% Q: State weights

% R: Penalization/weights on rate of change in u, uk - uk-1

% Ru: Control input u penalization/weights

%

% Output:

% J: objective function cost

%

%% Nonlinear MPC design parameters

%% Cost Calculation

% Set initial plant states, controller output and cost

xk = x;

uk = u(1);

J = 0;

% Loop through each prediction step

for ct=1:Np

% Obtain plant state at next prediction step

% xk1 = GTStateFcnDT(xk,uk);
xk1 = state_sim(xk,uk);

% yk1 = GTOutputFcn(xk1,uk);

% Accumulate state tracking cost from x(k+1) to x(k+Np)

% J = J + (yk1-xref)'*Q*(yk1-xref);
J = J + (xk1-xref)'*Q*(xk1-xref);

% Accumulate rate of change cost from u(k) to u(k+Np-1)

if ct==1

J = J + (uk-u0)'*R*(uk-u0) + uk'*Ru*uk;

else

J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + uk'*Ru*uk;

end

% Update xk and uk for the next prediction step

xk = xk1;

if ct<Np

uk = u(ct+1);

end

end