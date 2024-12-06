%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 5th order upwind SBP operator           %%%
%%%                                         %%%
%%% Derived by Ken Mattsson                 %%%
%%%                                         %%%
%%% Set number of unknowns m, grid step h   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
H=diag(ones(m,1),0);
H(1:4,1:4)=[0.251e3 / 0.720e3 0 0 0; 0 0.299e3 / 0.240e3 0 0; 0 0 0.211e3 / 0.240e3 0; 0 0 0 0.739e3 / 0.720e3;];
H(m-3:m,m-3:m)=fliplr(flipud(H(1:4,1:4)));
H=H*h;
HI=inv(H);
   
Qp=(1/20*diag(ones(m-2,1),-2)-1/2*diag(ones(m-1,1),-1)-1/3*diag(ones(m,1),0)+1*diag(ones(m-1,1),+1)-1/4*diag(ones(m-2,1),+2)+1/30*diag(ones(m-3,1),+3));

Q_U = [-0.1e1 / 0.120e3 0.941e3 / 0.1440e4 -0.47e2 / 0.360e3 -0.7e1 / 0.480e3; -0.869e3 / 0.1440e4 -0.11e2 / 0.120e3 0.25e2 / 0.32e2 -0.43e2 / 0.360e3; 0.29e2 / 0.360e3 -0.17e2 / 0.32e2 -0.29e2 / 0.120e3 0.1309e4 / 0.1440e4; 0.1e1 / 0.32e2 -0.11e2 / 0.360e3 -0.661e3 / 0.1440e4 -0.13e2 / 0.40e2;];

Qp(1:4,1:4)=Q_U;
Qp(m-3:m,m-3:m)=flipud( fliplr(Q_U(1:4,1:4) ) )'; %%% This is different from standard SBP

Qm=-Qp';

e_1=zeros(m,1);e_1(1)=1;
e_m=zeros(m,1);e_m(m)=1;

Dp=HI*(Qp-1/2*e_1*e_1'+1/2*e_m*e_m') ;

Dm=HI*(Qm-1/2*e_1*e_1'+1/2*e_m*e_m') ;



