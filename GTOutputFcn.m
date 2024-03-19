function y = GTOutputFcn(x,u)

 

global N tauGT TauL m phi

persistent Nv Tv tauL

 

D = 12;

% if x(2) ~= 0

 

    Ktau = tauGT/(x(2)*mean(m)*(D-6)*1e-3/0.08);

   

    tauGTsub = x(2)*(D-6)*Ktau*mean(m)*1e-3/0.08;

 

    dtau = tauGTsub - TauL;

 

    if 242 < dtau < 312 && dtau >~ 0 && dtau < 312 && dtau > 242

        Nv = 3200;

        Tv = tauGT;

        y = Nv;

        tauL = TauL;

    else

        pctT = dtau*100/Tv;

        if pctT >= 100

            dpct = pctT - 100;

            pctN = pctT + (dpct/100);

            y = Nv*pctN/100;

        elseif dtau < 0

            pctN = 100 + pctT;

            y = Nv*pctN/100;

        elseif 0 < pctT < 100 && pctT > 0 && pctT < 100 && phi < 0.55 || TauL > tauL

            pctN =  100 - (pctT/100);

            y = Nv*pctN/100;

        elseif 0 < pctT < 100 && pctT > 0 && pctT < 100 && phi > 0.55 || TauL < tauL

            pctN =  100 + (pctT/100);

            y = Nv*pctN/100;

        end

   

    end

% else

%     y = N;

% end

 

end