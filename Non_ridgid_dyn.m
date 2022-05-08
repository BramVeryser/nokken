clear; close all; 

% Constanten
cam = load('Geometrie_e.mat');
%rise pas aan
rise=3;
k_f = 10; 
k_s = 0;
m = 3;%18?
dzeta = 0.1;
omega = pi;
omega_n = sqrt( (k_f+k_s)/m);
% welke rise?
if (rise == 1) % waardes uit tabel slide 24
    beta = 75; % aanpassen voor verschillende rise
    N = 4;
    Q = 24*35;
elseif (rise == 2)
    beta = 60;
    N = 3;
    Q = 60;
elseif (rise == 3)
    beta = 80;
    N = 3; 
    Q = 60;
end
t_n = 2*pi/omega_n;
t_1 = 2*pi/180 * beta / omega;
lambda = t_1/t_n; % lambda is zeer klein <<10 dus benadering niet geldig
dzeta*lambda - 0.75 %moet groter zijn dan nul anders is benadering slecht
% het ziet er naar uit dat we numeriek gaan moeten uitwerken 
% k_f kiezen zodat er aan de benaderingsvoorwaarde voldaan wordt
k_f = 2*pi*m*.75/t_1/dzeta -k_s


theta = cam.S;
theta_dot  = cam.V;
t= 0:0.001:t_1; % in seconden
tau =  0:(1/800):10;%t/t_1;  %

A_max = Q/(2*pi*lambda)^N


%% numerieke analyse
teller = (2*pi*lambda)^2;
noemer = [1, 2*dzeta*(2*pi*lambda), (2*pi*lambda)^2];
sys = tf(teller, noemer);

theta_3 = theta(200e2:280e2);
theta_init = [theta(200e2) theta_dot(200e2)];
lsim(sys,theta_3/30-1,tau, theta_init)



    