function y = RK4(x, k, M)
    k1 = M*x;
    k2 = M*(x + (k./2).*k1);
    k3 = M*(x + (k./2).*k2);
    k4 = M*(x + k.*k3);
    y = x + (k./6).*(k1 + 2.*k2 + 2.*k3 + k4);