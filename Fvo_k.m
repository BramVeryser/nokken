clear; close all; clc

out = load('Geometrie.mat');

X = out.S;
V = out.V;
A = out.A;
alpha = out.pressure_angle;
omega = out.w;
F_func = out.extload;
m = out.mass;
rho = out.roc_pitch;

%plot(rho)
%min(abs(rho))
N_tot = (F_func + m*(omega^2)*(A*10^-3))./cos(alpha);

%[N_min,I_Nmin] = min(N_tot);

I_Nmin = 24000;
N_min = N_tot(I_Nmin);
Fv0 = zeros(1,201);
N_max = zeros(1,201);
Norm = zeros(1,201);
for k = 0:0.1:20
    Fv0(round(k*10+1)) = -N_min*cos(alpha(I_Nmin)) - k*X(I_Nmin);
    N_max(round(k*10+1)) = max((F_func + Fv0(round(k*10+1)) + X*k + m*(omega^2)*(A*10^-3))./cos(alpha));
    N_nieuw = (F_func + Fv0(round(k*10+1)) + X*k + m*(omega^2)*(A*10^-3))./cos(alpha);
    Norm(round(k*10+1)) = norm(N_nieuw,1)/36000;
end
figure()
plot(Fv0)

figure()
plot(N_max)

figure()
plot(Norm)




