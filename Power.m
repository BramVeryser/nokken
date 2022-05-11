clear; close all; clc

out = load('Geometrie_e0.mat');

X = out.S*10^-3;
V = out.V;
A = out.A;
alpha = out.pressure_angle;
omega = out.w;
F_func = out.extload;
m = out.mass;
rho = out.roc_pitch;
N_tot = out.normalforce_tot;
e = out.exc*10^-3;
bcr = out.bcr*10^-3;
rof = out.rof*10^-3;

R0 = bcr + rof;
R = X + R0;
P_e0 = N_tot.*sin(alpha).*R.*omega;
figure()
plot(P_e0)


test = load('Geometrie_e.mat');

Xe = test.S*10^-3;
V = out.V;
A = out.A;
alphae = test.pressure_angle;
omega = out.w;
F_func = out.extload;
m = out.mass;
rho = out.roc_pitch;
N_tote = test.normalforce_tot;
e = test.exc*10^-3;
bcre = test.bcr*10^-3;
rofe = test.rof*10^-3;

R0e = bcre + rofe;
R_no_exc = sqrt(test.xpitch.^2 + test.ypitch.^2);%*10^(-3);

P_e = N_tote.*(e*cos(alphae)+sin(alphae).*(sqrt(R0e^2-e^2)+Xe))*omega;
P_av2 = mean(P_e);
figure()
plot(P_e)


M = P_e/omega;
M_av = P_av2 ./ test.w;
M_delta = M_av - M;
A = cumtrapz(test.theta, M_delta);

figure()
plot(M)
figure()
plot(A)