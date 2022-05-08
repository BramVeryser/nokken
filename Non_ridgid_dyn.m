clear; close all; 

% Constanten
cam = load('Geometrie_e.mat');
%rise pas aan
rise=1;
k_f = 10; 
k_s = 10;
m = 1;%18?
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

theta = cam.S;

t= 0:0.001:t_1; % in seconden
tau =  0:0.001:1;%t/t_1;  %

A_max = Q/(2*pi*lambda)^N


    