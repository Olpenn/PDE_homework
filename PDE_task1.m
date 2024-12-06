rho = 1;
c = 1;
beta = 0;
m = 51;
z = zeros(m);
h = 2/(m-1);
al = rho.*c;
ar = -rho.*c;

%%
C1 = diag(repelem(rho.*(c.^2), m));
C2 = diag(repelem(1./rho, m));
C_inv = [C1, z; z, C2];
D = [diag(repelem(beta, m)), z; z, z];


[Dp, Dm, HI] = upwind7(h, m);
Hinv = [HI, z; z, HI];
Dx = [z, Dp; Dm, z];

% dir

L = zeros(2, 2*m);
L(1, m+1) = 1;
L(1, end) = 1;

P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
Mdir = -P*C_inv*(Dx + D)*P;

% char

L = zeros(2, 2*m);
L(1, 1) = 1;
L(2, m) = 1;
L(1, m+1) = al;
L(2, end) = ar;

P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
Mchar = -P*C_inv*(Dx + D)*P;

lambda_dir = eig(Mdir);
lambda_char = eig(Mchar);
figure
plot(real(lambda_dir.*h));
figure
plot(real(lambda_char.*h));

rk_stab = @(x, lambda) 1 + lambda.*x + ((lambda.*x).^2)/2 + ((lambda.*x).^3)/6 + ((lambda.*x).^4)/24;

[lmax_dir, id] = max(abs(real(lambda_dir)));
[lmax_char, ic] = max(abs(real(lambda_char)));

CFL_arr = 0:0.01:2;
figure
plot(rk_stab(CFL_arr.*h, real(lambda_dir(id))), color='blue');
figure
plot(rk_stab(CFL_arr.*h, real(lambda_char(ic))), Color='red');
cfl = min(rk_stab(CFL_arr.*h,lmax_char));
%%

theta1 = @(x,t) exp(-(((x-t)./0.1).^2));
theta2 = @(x,t) -exp(-(((x+t)./0.1).^2));

tspan = [0, 1.8];
xspan = [-1, 1];
m_arr = [101];
err7 = zeros(5, 1);
conv7 = zeros(4, 1);
err6 = zeros(5, 1);
conv6 = zeros(4, 1);

for i = 1:size(m_arr, 2)
    m = m_arr(i);
    h = 2/(m-1);
    z = zeros(m);
    xtemp = xspan(1):h:xspan(2);
    x = [theta1(xtemp, 0) - theta2(xtemp, 0), theta1(xtemp, 0) + theta2(xtemp, 0)]';
    x_ex = [theta2(xtemp, 2-tspan(2)) - theta1(xtemp, 2-tspan(2)), theta1(xtemp, 2-tspan(2)) + theta2(xtemp, 2-tspan(2))]';
    cfl = 0.05;
    
    % 7

    C1 = diag(repelem(rho.*(c.^2), m));
    C2 = diag(repelem(1./rho, m));
    C_inv = [C1, z; z, C2];
    D = [diag(repelem(beta, m)), z; z, z];

    [Dp, Dm, HI] = upwind7(h, m);
    Hinv = [HI, z; z, HI];
    Dx = [z, Dp; Dm, z];
    L = zeros(2, 2*m);
    L(1, 1) = 1;
    L(2, m) = 1;

    P = eye(2*m) - Hinv*L'*inv(L*Hinv*L')*L;
    M7 = -P*C_inv*(Dx + D)*P;
    tcurr = 0;
    while tcurr <= tspan(2)
        x_new = RK4(x, cfl.*h, M7);
        tcurr = tcurr + cfl.*h;
        x_new(1) = 0;
        x_new(m) = 0;
        x = x_new;
    end
    figure
    plot(xtemp, x(1:m), 'red', xtemp, x(m+1:end))
    figure
    plot(xtemp, x_ex(1:m), 'red', xtemp, x_ex(m+1:end))
    error7 = sqrt(h).*norm(x_ex - x);
    if i > 1
        q7 = log10(err7(i-1)./error7)./log10(m/m_arr(i-1));
        conv7(i-1) = q7;
    end
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
        x_new(m+1) = 0;
        x_new(end) = 0;
        x = x_new;
    end
    error6 = sqrt(h).*norm(x_ex - x);
    if i > 1
        q6 = log10(err7(i-1)./error6)./log10(m/m_arr(i-1));
        conv6(i-1) = q6;
    end
end


