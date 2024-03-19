close all;clc;
figure
subplot(1,2,1)
plot(tspan,rHistory,'*-',tspan,yHistory,'LineWidth',2)
xlabel("Time")
ylabel("RPM")
legend(["RPM\_ref" "RPM\_GT"],'Location','southeast')
subplot(1,2,2)
plot(tspan,Mf,tspan,MdotfHistory,'c*','LineWidth',2)
xlabel("Time")
ylabel("Mass\_Flow")
legend(["Mass\_Flow\_NMPC" "Mass\_Flow\_GT"],'Location','southeast')
