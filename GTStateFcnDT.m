function xk1 = GTStateFcnDT(xk,uk)

Ts = 0.04;

xk1 = xk(:);

nu = 10;  % Number of integration time steps for Euler method

dt = Ts/nu;

uk1 = uk(:);

for i = 1:nu
    xk1 = xk1 + dt*GTStateFcnCT(xk1,uk1);
end

end