function [QF,Fsr,Fsrt] = SimFSR(Tref,FF_FSR,FF_FSRT,Cv,Sg,dp,Anm)

 

global mdotF T0x N

 

% initial conditions

Kp1 = 4.2;

Ki1 = 0.001;

Kp2 = 0.1;

Ki2 = 0.0001;

Kd = 10.5;

 

ff0 = 0;

ff0T = 0;

TNR = 104;

 

T = 1500;

dt = 0.08;

tspan = 0:dt:T;

 

for k=1:numel(tspan)

 

    % FSRN

    TNH = N*100/3200;

   

    if (TNR-TNH) ~= 4 && (TNR-TNH) < 20 && (TNR-TNH) > 0

        FSRN = (TNR-TNH)*12 + 16;

    else

        FSRN = 70;

    end

 

    % FSRT

    FSRT_PS = FF_FSRT*100/70;

    TTXM = T0x - Tref + Tref*FSRT_PS/100;

   

    % FF PID FSRT

    [~, xffT] = ode45(@(t,x) PID(t,x,Tref,TTXM,Kp2,Ki2,Kd), [k*0.0001 (k+1)*0.0001], ff0T);

    PIRegT = xffT(end,1);

    ff0T = xffT(end,:);

    FF_FSRT = FF_FSRT + PIRegT;

 

    % Aux FSR

    FSR = min([FF_FSRT,FSRN]);

 

    % FSR

    Arfsr = (((mdotF/Cv)^2)*Sg)/dp;

    FSRr = Arfsr*70/Anm;

 

    % FF PI FSR

    [~, xff] = ode45(@(t,x) PI(t,x,FSRr,FF_FSR,Kp1,Ki1), [k*0.0001 (k+1)*0.0001], ff0);

    PIReg = xff(end,1);

    ff0 = xff(end,:);

    FF_FSR = FSR + PIReg;

 

    % calc fsr

    Afsr = Anm*FF_FSR/70;

    Qf = Cv*sqrt(Afsr*dp/Sg);

 

end

 

QF = Qf;

Fsr = FF_FSR;

Fsrt = FF_FSRT;

end