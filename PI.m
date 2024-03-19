function dxdt = PI(t, x, FSRr, FSRfb, Kp, Ki)

    dxdt = zeros(1,1);

    % construction of PID

    r = FSRr;                           % reference signal

    e = FSRfb - r;                      % error signal

    u = - Kp*e - Ki*x(1);   % the PID thing

    dxdt(1) = e*x(1) + u;               % for integral action in PID

end