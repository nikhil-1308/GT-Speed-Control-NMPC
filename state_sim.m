function x1 = state_sim(x,u)
x1 = x;
T = 200;
dt = 0.08;
tspan = 0:dt:T;

for i=1:numel(tspan)

    x1 = GTStateFcnDT(x1,u);
    
end