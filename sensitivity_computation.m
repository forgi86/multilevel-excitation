%% Computation of the sensitivity equations which are integrated together with the nonlinear model in order to obtain the information matrix F %%

addpath('util');
clear;

%% SENSITIVITY COMPUTATIONS %%
syms F VR Ca Ca0 k0 E R TR T0 lambda FJ VJ TCin TJ UAJ rho cp VJ rhoJ cJ;
x = [Ca;TR;TJ];

theta = [k0 E lambda UAJ];

f1 = F/VR*(Ca0 - Ca) - Ca*k0*exp(-E/(R*TR));
f2 = F/VR*(T0 - TR) - lambda*Ca*k0*exp(-E/(R*TR))/(rho*cp) -(UAJ/(VR*rho*cp))*(TR-TJ);
f3 = FJ/(VJ)*(TCin - TJ) + (UAJ/(VJ*rhoJ*cJ))*(TR - TJ);
f = [f1;f2;f3];

g = [TR];


Jfx = simplify(jacobian(f,x));
Jp = simplify(jacobian(f,theta));

n = length(x);
m = length(theta);
Wx = symMat([n m],'w');

Wxdot = simple(Jfx*Wx + Jp);

Wy = jacobian(g,x)*Wx;