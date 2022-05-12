clear; close all; 

% Constanten
cam = load('Geometrie_e.mat');
%rise pas aan
rise=2;
%k_f = 10; moeten we bepalen mogen we niet vanuit gaan
k_s = 0;
m = 18;
dzeta = 0.1;
omega = pi;
%omega_n = sqrt( (k_f+k_s)/m);

% welke rise? -> 2 om kf te bepalen
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
%% Eerst kf bepalen zodat benaderende oplossing juist is. 
%Zie hfdst 9 slide 30: dzeta*lambda > 0.75  

lamda_min = 0.75/dzeta; % lambda > 7.5

%Slide 11: lamda is t1/t_n, t1 is de tijd waarin de rise/fall gebeurt
%t_n = 2pi/omega_n
%dus lamda = t1*(omega_n/(2pi))( > 5.53)
%en omega_n = sqrt( (k_f+k_s)/m) met ks = 0
%ook bij de rise/fall in de korste tijd moet de voorwaarde op lamda voldaan zijn
%(Bij de andere is het dan ook voldaan omdat de rest constant is)
%Dus hier is het de tweede rise 
t_1 = 2*pi/180 * beta / omega;

t_1_A = 2*pi/180 * 75 / omega;

t_1_C = 2*pi/180 * 80 / omega;

k_f_exact = m*4*pi^2*lamda_min^2/(t_1^2);   %k_f moet groter zijn dat dit getal ([N/mm] denk ik lijkt mij het meest logisch)
k_f = 90000;                      % neem 90000

t_n = 2*pi/(sqrt(k_f/m));

%check of de rest klopt moeten allemaal groter zijn dan 0.75
assert(t_1/t_n*dzeta > 0.75,"k_f is te klein")
assert(t_1_A/t_n*dzeta > 0.75,"k_f is te klein")
assert(t_1_C/t_n*dzeta > 0.75,"k_f is te klein")

%% Numerieke analyse 
%kies opnieuw rise 2 (de enige voorwaarde hier is dat die door een dwell moet gevolgd worden, 3 zou dus ook gaan)
lambda = t_1/t_n;
teller = (2*pi*lambda)^2;
noemer = [1,2*dzeta*(2*pi*lambda),(2*pi*lambda)^2];
sys = tf(teller,noemer);

%simulatie
Ts = 0.001;
tau = (0:Ts:8)/6;  %tau is 1 als heffing gedaan is maar we simuleren nog iets verder
theta = cam.S(1,12000:20000)/15-1; %de -1 is om er voor te zorgen dat we starten bij theta = 0 en eindigen bij theta = 1
theta0 = 0;
theta_dot0 = 0;
[A,B,C,D] = tf2ss(teller,noemer);
X0=[1/C(2)*theta_dot0;1/C(2)*theta0];
lsim(A,B,C,D,theta,tau,X0)
gamma = lsim(A,B,C,D,theta,tau,X0);
gamma_dot = diff(gamma')./Ts;
gamma_dot_dot = diff(gamma',2)./Ts^2;
gamma_dot_dot = [gamma_dot_dot,gamma_dot_dot(end),gamma_dot_dot(end)];

verschil_single=gamma' - theta;

%figuur veschil
figure()
plot(tau(6000:8000),verschil_single(6000:8000))

%omhullende bepalen A1 via slide 13
x0 = gamma(6000)-1;
lambda_d = lambda*sqrt(1-dzeta^2);
v0 = (gamma(6001)-gamma(5999))/(tau(6001)-tau(5999));  %snelheid door numeriek afleiden
A1 = sqrt(((x0*2*pi*lambda_d)^2+(v0+dzeta*2*pi*lambda*x0)^2)/((2*pi*lambda_d)^2));

omhullende_1 = A1*exp(-dzeta*2*pi*lambda*(tau(6000:8000)-1));
omhullende_2 = -A1*exp(-dzeta*2*pi*lambda*(tau(6000:8000)-1));

figure()
plot(tau(6000:8000), verschil_single(6000:8000), 'b', tau(6000:8000), omhullende_1, 'b', tau(6000:8000), omhullende_2, 'b')

%% Benadering
%ik heb enkel A1_tilde gedaan en ni de volledige functie aangezien ze
A1_tilde = (Q/(2*pi*lambda)^N)*sqrt(1/(1-dzeta^2));

epsilon = (A1-A1_tilde)/A1;
% M=2*pi*lambda;
% b_2=-(dzeta*Q)/M;
% b_1=(4*dzeta^2*Q-Q)/M^2;
% b_0=1/M^2 + 2*Q*dzeta/M^3-(2*dzeta*(4*dzeta^2*Q-Q))/M^3;
% %b_0=2*Q*dzeta/M^3-(2*dzeta*(4*dzeta^2*Q-Q))/M^3;
% x0 = (Q*(tau-1).^N)/factorial(N) + b_0 + b_1*(tau-1) + b_2*(tau-1).^2-1;
% v0 = Q*(tau-1)^2/2+b_1+2*b_2*(tau-1);
% theta_benadering = (Q*(tau-1)^N)/factorial(N)+1;
% A1 = sqrt(((x0*2*pi*lambda_d)^2+(v0+dzeta*2*pi*lambda*x0)^2)/((2*pi*lambda_d)^2))

%% Invloed op kracht
out = load('Geometrie_e.mat');

X = out.S;
V = out.V;
A = out.A;
alpha = out.pressure_angle;
omega = out.w;
F_func = out.extload;
m = out.mass;
rho = out.roc_pitch;
k_v = out.springconstant;
F_v = out.springpreload+10;


%plot(rho)
%min(abs(rho))
N_tot = (F_func + X*k_v +F_v+ m*(omega^2)*(A*10^-3))./cos(alpha);
%N_tot_dynamica = N_tot(12000:20000)+(k_v*15*verschil_single)./cos(alpha(12000:20000));
N_tot_dynamica = (F_func(12000:20000) + gamma'*15*k_v + F_v + m*(omega^2)*gamma_dot_dot*15)./cos(alpha(12000:20000));
figure()
plot(N_tot(12000:20000))
hold on
plot(N_tot_dynamica)

figure()
plot(N_tot_dynamica-N_tot(12000:20000))

min(N_tot_dynamica)

%% multirise
T = 2*pi/cam.w;
lift = cam.S*0.001/0.03;
grid on
lift = repmat(lift,1,25);

tau = (0:1:899999)*25/900000;

lambda = T/t_n;
teller = (2*pi*lambda)^2;
noemer = [1,2*dzeta*(2*pi*lambda),(2*pi*lambda)^2];
sys = tf(teller,noemer);

gamma_multi = lsim(sys,lift,tau);
gamma_multi_ss = gamma_multi(864001:900000)';
gamma_dot_dot_ss = diff(gamma_multi_ss,2)/Ts^2;
gamma_dot_dot_ss = [gamma_dot_dot_ss,gamma_dot_dot_ss(end),gamma_dot_dot_ss(end)];
%Volledige cyclus
plot(tau(864001:900000),gamma_multi(864001:900000),"linewidth",1.7)
grid on
title("Full 20'th cycle \gamma")
xlabel("\tau")
ylabel("\gamma")
%trillingen
verschil_multi = gamma_multi'-lift;
figure()
plot(tau(864001:900000),verschil_multi(864001:900000),"linewidth",1.7)
title("Difference motion law and \gamma_{multi}")
grid on
xlabel("\tau")
ylabel("\gamma - S")
figure()
plot(verschil_single,"linewidth",1.7)
hold on
plot(verschil_multi(876000:884000),"linewidth",1.7)
legend("\gamma_{single} - S","\gamma_{multi} - S")
title("Difference \gamma_{single} and \gamma_{multi} on critical rise")
grid on
xlabel("\tau")
ylabel("\Delta\gamma")
figure()
plot(gamma,"linewidth",1.7)
hold on
plot(gamma_multi(876000:884000),"linewidth",1.7)
title("gamma/gammamulti")
grid on
%% Invloed kracht multi
out = load('Geometrie_e.mat');

X = out.S;
V = out.V;
A = out.A;
alpha = out.pressure_angle;
omega = out.w;
F_func = out.extload;
m = out.mass;
rho = out.roc_pitch;
k_v = out.springconstant+11;
F_v = out.springpreload+10+110;


%plot(rho)
%min(abs(rho))
N_tot = (F_func + X*k_v +F_v+ m*(omega^2)*(A*10^-3))./cos(alpha);
%N_tot_dynamica = N_tot(12000:20000)+(k_v*15*verschil_single)./cos(alpha(12000:20000));
N_tot_dynamica = (F_func + gamma_multi_ss*30*k_v + F_v + m*(omega^2)*gamma_dot_dot_ss*30)./cos(alpha);
figure()
plot(N_tot,"linewidth",1.7)
xlabel("cam angle [deg]")
ylabel("N_{tot} [N]")
title("Total normal force")
hold on
plot(N_tot_dynamica,"linewidth",1.7)
legend("N(S)","N(\gamma)")
% 
% figure()
% plot(N_tot_dynamica-N_tot)
% figure()
% plot(gamma_dot_dot_ss, "linewidth",1.7)
% title("Acceleration \gamma")
% grid on
% 
% xlabel("\tau")
% ylabel("d^2\gamma/d\tau^2")
% min(N_tot_dynamica)

% figure()
% plot(tau(12001:18000),gamma(1:6000))
% plot(tau(12001:18000),gamma_multi(12001:18000))
% 
% gamma_corrected = (gamma_multi(18001+900000-36000:20000+900000-36000)-0.5)*2;
% % plot(tau(12001:18000),gamma_corrected)
% 
% %singlerise en multrise vergelijken
% figure()
% plot(tau(18001:20000),gamma(6001:8000)-lift(18001:20000))
% hold on
% plot(tau(18001:20000),gamma_corrected-lift(18001:20000))
% 
% figure()
% plot(tau(18001:20000),gamma(6001:8000)-gamma_corrected)


%% Oude code
% lambda = t_1/t_n; % lambda is zeer klein <<10 dus benadering niet geldig
% dzeta*lambda - 0.75 %moet groter zijn dan nul anders is benadering slecht
% % het ziet er naar uit dat we numeriek gaan moeten uitwerken 
% % k_f kiezen zodat er aan de benaderingsvoorwaarde voldaan wordt
% k_f = 2*pi*m*.75/t_1/dzeta -k_s
% 
% 
% theta = cam.S;
% theta_dot  = cam.V;
% t= 0:0.001:t_1+0.8889; % in seconden
% tau =  0:0.001:16;%t/t_1;  %
% 
% A_max = Q/(2*pi*lambda)^N
% 
% 
% %% numerieke analyse
% teller = (2*pi*lambda)^2;
% noemer = [1, 2*dzeta*(2*pi*lambda), (2*pi*lambda)^2];
% sys = tf(teller, noemer);
% 
% theta_3 = theta(200e2:360e2);
% theta_init = [theta(200e2) theta_dot(200e2)];
% lsim(sys,theta_3/30-1,tau, theta_init)
% y = lsim(sys,theta_3/30-1,tau, theta_init);
% y_dot = diff(y');
% y_dot_dot = diff(y',2);
% figure()
% plot(y_dot)
% title("y_{dot}")
% figure()
% plot(y_dot_dot)
% title("y_{dotdot}")
% 


    