function dxdt = PID(t, x, FSRTr, FSRTfb, Kp, Ki, Kd)

    dxdt = zeros(1,1);

    % construction of PID

    r = FSRTr;                           % reference signal

    e = FSRTfb - r;                      % error signal

    u = - Kp*e - Ki*x(1)- Kd*x(1);   % the PID thing

    dxdt(1) = e*x(1) + u;               % for integral action in PID

end

 