function dxdt = GTStateFcnCT(x,u)

 

global tauGT N TauL Ain alpha Cf m phi mdotF T0x

 

Parameters;

 

% inputs

phi = u;

Patm = P0;

Tatm = T0;

 

U = (pi*N*D)/60;

%% Compressor Stage 

for i=1:numel(len)-1

    Lc = L - dx;

    y(i) = tan_theta*Lc;

    if i==1

        Vin = (D - 6 + 2*y(i))^2*(pi/4)*dx;

    end

    Vx(i) = (D - 6 + 2*y(i))^2*(pi/4)*dx;

    P2(i) = ((V0/Vx(i))^gam)*P0;

    T2(i) = ((P2(i)/P0)^(1-(1/gam)))*T0;

    roh2(i) = roh0 + ((R*T2(i))/(P2(i)*Vx(i)))*1e-4;

 

    % velocity equations

    Nu2(i) = sqrt(2*(P0 + 0.5*roh0*Nu0^2 - 4*P2(i))/roh2(i));

    Nuf(i) = (U - cos(bet))*Vr*tan(alpha);

    NuT(i) = Nuf(i) + Nu2(i)*0.08;

   

    % update step

    P0 = P2(i);

    T0 = T2(i);

    roh0 = roh2(i);

    V0 = Vx(i);

    Nu0 = Nu2(i);

    L = Lc;

   
end

%% Combustion Stage 

% parameters

CmbRt = phi*0.1e3;

Vc1 = Vx(end);

Vc2 = Vc1 + 3;

dt = 0.08;

T = 1;

tspan = 0:dt:T;

Cv = 10.6;

 

% initial conditions

Roh0 = roh2(end);

Tc0 = T2(end);

Pc0 = P2(end);

Qp = CmbRt*Pc0;

NuC0 = NuT(end);

 

for k=1:numel(tspan)

 

    Roh2(k) = Roh0 + (Roh0/(1 + Qp));

    Tc2(k) = Tc0 - 0.08*Qp;

    dS(k) = Cv*log(Tc2(k)/Tc0) + R*log(Vc2/Vc1);

 

    % final temp

    Tcs(k) = Qp/dS(k);

    if k>=2

        Tc2(k) = Tc2(k) + Tcs(k);

    end

    Tct(k) = Tc0 + Tc2(k)*dt;

    NuC(k) = NuC0 + sqrt(3*1.38e-2*Tct(k)/Roh2(k));

 

    % update step

    Roh0 = Roh2(k);

    Tc0 = Tc2(k) + 0.08*Qp;

    NuC0 = NuC(k);

 

end

%% Turbine Stage 

mdotF = Qp/(Cf*Eff);

mdotA = mdotF/phi;

Af = mdotA*R*Tatm/(Patm*Vin);

alpha = atand(Af/Ain);

if alpha >= 11

    alpha = 10;

    Ain = 14.67;

    Cf = (tand(alpha)*Ain*Patm*Vin*Eff)*1e5/(R*Tatm*phi*Qp);

 

    % post calculation

    mdotF = Qp/(Cf*Eff);

    mdotA = mdotF/phi;

    Af = mdotA*R*Tatm/(Patm*Vin);

    alpha = atand(Af/Ain);

    Ain = Af;

elseif alpha < 11

    Ain = Af;

end

 

clear i

 

% parameters

Lx = L;

V0x = Vc2;

P0x = Pc0;

T0x = Tct(end);

ron0x = Roh2(end);

Nu0x = NuC(end);

DL = 1.5;

AL = (pi/4)*DL^2;

FL = TauL/DL^2;

M = 100;

m = 1:M/numel(len):M;

 

% numel(m)

for i=1:numel(len)-1

    Lcx = Lx + dx;

    y(i) = tan_theta*Lcx;

    Vxx(i) = (D - 6 + 2*y(i))^2*(pi/4)*dx;

    P2x(i) = ((Vc2/Vxx(i))^gam)*P0x;

    T2x(i) = ((P2x(i)/P0x)^(1-(1/gam)))*T0x;

    roh2x(i) = ron0x - ((R*T2x(i))/(P2x(i)*Vxx(i)))*1e-5;

 

    % velocity equations

    Nu2x(i) = sqrt(2*(P0x + 0.5*ron0x*Nu0x^2 - 2*P2x(i))/roh2x(i));

   

    NuTx(i) =  Nu2x(i);

 

    Fp(i) = (1/3)*roh2(i)*NuTx(i)^2;

    I(i) = 0.5*m(i)*(D - 6 + 2*y(i))^2;

    F(i) = Fp(i) - (I(i) + FL)*1e-2 + Kx;

    P(i) = P2x(i) - 1e-3*(1e-3*Fp(i) - (1e-3*I(i)/((D - 6 + 2*y(i))^2*(pi/4)) + FL/AL));

 

    % update step

    P0x = P2x(i);

    T0x = T2x(i);

    ron0x = roh2x(i);

    Vc2 = Vxx(i);

    Nu0x = Nu2x(i);

    Lx = Lcx;

   

end

%% Dynamics 

% sys = tf(mean(F),[m(end) (I + FL) 10]); %step(sys)

[A,B,~,~] = tf2ss(1,[mean(m) mean(I)+FL Kx]);

 

dxdt = A*x + B*mean(Fp);

 

tauGT = mean(F)*(D-6)*1e-3;

 

end