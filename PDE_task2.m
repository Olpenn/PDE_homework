rho1 = 1;
c1 = 1;
rho2 = 2;
c2 = 2;
beta = 0;
al = rho1.*c1;
ar = -rho2.*c2;

theta1 = @(x,t) exp(-(((x-t)./0.1).^2));
theta2 = @(x,t) -exp(-(((x+t)./0.1).^2));

tspan = [0, 2.3];
xspan = [-2, 4];
x_int = 1;
m_arr = [201, 401];
err7 = zeros(5, 1);
conv7 = zeros(4, 1);
err6 = zeros(5, 1);
conv6 = zeros(4, 1);

for i = 1:size(m_arr, 2)
    m = m_arr(i);
    h = (xspan(2) - xspan(1))/(m-1);
    z = zeros(m);
    xtemp = xspan(1):h:xspan(2);
    x = [theta1(xtemp, 0) - theta2(xtemp, 0), theta1(xtemp, 0) + theta2(xtemp, 0)]';
    cfl = 0.05;
    
    % 7
    i_int = round(((x_int - xspan(1))./h) + 1);
    C1_vec = repelem(rho1.*(c1.^2), m);
    C1_vec(i_int:end) = repelem(rho2.*(c2.^2), size(C1_vec(i_int:end), 2));
    C1 = diag(C1_vec);
    C2_vec = repelem(1./rho1, m);
    C2_vec(i_int:end) = repelem(1./rho2, size(C2_vec(i_int:end), 2));
    C2 = diag(C2_vec);
    C_inv = [C1, z; z, C2];
    D = [diag(repelem(beta, m)), z; z, z];

    [Dp, Dm, HI] = upwind7(h, m);
    Hinv = [HI, z; z, HI];
    Dx = [z, Dp; Dm, z];
    L = zeros(2, 2*m);
    L(1, 1) = 1;
    L(2, m) = 1;
    L(1, m+1) = al;
    L(2, end) = ar;

    P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
    M7 = -P*C_inv*(Dx + D)*P;
    tcurr = 0;
    while tcurr <= tspan(2)
        x_new = RK4(x, cfl.*h, M7);
        tcurr = tcurr + cfl.*h;
        x_new(1) = -al.*x_new(m+1);
        x_new(m) = -ar.*x_new(end);
        x = x_new;
    end
    figure
    plot(xtemp, x(1:m), 'red', xtemp, x(m+1:end))
    % 6

    [D1, HI] = cent6(h, m);
    Hinv = [HI, z; z, HI];
    Dx = [z, D1; D1, z];

    P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
    M7 = -P*C_inv*(Dx + D)*P;
    tcurr = 0;
    while tcurr <= tspan(2)
        x_new = RK4(x, cfl.*h, M7);
        tcurr = tcurr + cfl.*h;
        x_new(1) = -al.*x_new(m+1);
        x_new(m) = -ar.*x_new(end);
        x = x_new;
    end
    plot(xtemp, x(1:m), 'red', xtemp, x(m+1:end))
end
%%

rho1 = 1;
c1 = 1;
rho2 = 2;
c2 = 2;
beta1 = 0;
beta2 = 10;
al = rho1.*c1;
ar = -rho2.*c2;

theta1 = @(x,t) exp(-(((x-t)./0.1).^2));
theta2 = @(x,t) -exp(-(((x+t)./0.1).^2));

tspan = [0, 2.5];
xspan = [-4, 6];
x_int1 = 1;
x_int2 = 1.8;
x_comp = 0.2;
m_arr = [401];
err7 = zeros(5, 1);
conv7 = zeros(4, 1);
err6 = zeros(5, 1);
conv6 = zeros(4, 1);

for i = 1:size(m_arr, 2)
    m = m_arr(i);
    h = (xspan(2) - xspan(1))/(m-1);
    z = zeros(m);
    xtemp = xspan(1):h:xspan(2);
    x = [theta1(xtemp, 0) - theta2(xtemp, 0), theta1(xtemp, 0) + theta2(xtemp, 0)]';
    cfl = 0.05;
    
    % 7
    i_int1 = round(((x_int1 - xspan(1))./h) + 1);
    i_int2 = round(((x_int2 - xspan(1))./h) + 1);
    i_comp = round(((x_comp)./h) + 1);
    C1_vec = repelem(rho1.*(c1.^2), m);
    C1_vec(i_int1:(i_int1+i_comp)) = repelem(rho2.*(c2.^2), size(C1_vec(i_int1:(i_int1+i_comp)), 2));
    C1_vec(i_int2:(i_int2+i_comp)) = repelem(rho2.*(c2.^2), size(C1_vec(i_int2:(i_int2+i_comp)), 2));
    C1 = diag(C1_vec);
    C2_vec = repelem(1./rho1, m);
    C2_vec(i_int1:(i_int1+i_comp)) = repelem(1./rho2, size(C1_vec(i_int1:(i_int1+i_comp)), 2));
    C2_vec(i_int2:(i_int2+i_comp)) = repelem(1./rho2, size(C1_vec(i_int2:(i_int2+i_comp)), 2));
    C2 = diag(C2_vec);
    C_inv = [C1, z; z, C2];
    D_vec = repelem(beta1, m);
    D_vec(i_int1:(i_int1+i_comp)) = repelem(beta2, size(D_vec(i_int1:(i_int1+i_comp)), 2));
    D_vec(i_int2:(i_int2+i_comp)) = repelem(beta2, size(D_vec(i_int2:(i_int2+i_comp)), 2));
    D = [diag(D_vec), z; z, z];

    [Dp, Dm, HI] = upwind7(h, m);
    Hinv = [HI, z; z, HI];
    Dx = [z, Dp; Dm, z];
    L = zeros(2, 2*m);
    L(1, 1) = 1;
    L(2, m) = 1;
    L(1, m+1) = al;
    L(2, end) = ar;

    P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
    M7 = -P*C_inv*(Dx + D)*P;
    tcurr = 0;
    while tcurr <= tspan(2)
        x_new = RK4(x, cfl.*h, M7);
        tcurr = tcurr + cfl.*h;
        x_new(1) = -al.*x_new(m+1);
        x_new(m) = -ar.*x_new(end);
        x = x_new;
    end
    figure
    plot(xtemp, x(1:m), 'red', xtemp, x(m+1:end))
    % 6

    [D1, HI] = cent6(h, m);
    Hinv = [HI, z; z, HI];
    Dx = [z, D1; D1, z];

    P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
    M7 = -P*C_inv*(Dx + D)*P;
    tcurr = 0;
    while tcurr <= tspan(2)
        x_new = RK4(x, cfl.*h, M7);
        tcurr = tcurr + cfl.*h;
        x_new(1) = -al.*x_new(m+1);
        x_new(m) = -ar.*x_new(end);
        x = x_new;
    end
end