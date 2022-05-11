%this script is a replacement for the graphs on slide 38-39 of chapter 8
%this script will generate graphs, similar to the ones in the course text,
%that can be used for all cases
function [rho_min_beta] = gen_fig_Kloomok_Muffley_2(R0,beta_test,start_lift,end_lift,motionlaw)

%% input
% R0 = 60; %pitch radius in mm
% start_lift = 0; %start lift in mm
% end_lift = 15; %end lift in mm
% motionlaw = 4; 

%% calculations
beta_vec = 1:0.1:360;
rho_min = 0*beta_vec;

L0 = start_lift;
L1 = end_lift;
for i = 1:length(beta_vec)
    beta = beta_vec(i)*pi/180;
    theta = linspace(0, beta_vec(i), 100)*pi/180;
    
    x = theta/beta;
    
    if motionlaw==1 % dwell
            assert(L1-L0 == 0, 'The given input does not represent a dwell');
            S=L0*ones(1,size(x,2));
            V=0*ones(1,size(x,2));
            A=0*ones(1,size(x,2));

    elseif motionlaw==2 % minimal rms acceleration
            L=L1-L0;
            S=L*(3*x.^2-2*x.^3)+L0;
            V=L/beta*(6*x-6*x.^2);
            A=L/beta^2*(6-12*x);

    elseif motionlaw==3 % harmonische
            L=L1-L0;
            S=L*(1-cos(pi*x))/2+L0;
            V=L/beta*sin(pi*x)*pi/2;
            A=L/beta^2*cos(pi*x)*pi^2/2;

    elseif motionlaw==4 % volle cycloide
            L=L1-L0;
            S=L*(x-sin(2*pi*x)/2/pi)+L0;
            V=L/beta*(1-cos(2*pi*x));
            A=2*pi*L/beta^2*(sin(2*pi*x));

    elseif motionlaw==5 % 5th degree poly
            L=L1-L0;
            S = L0+L*(6*x.^5-15*x.^4+10*x.^3);
            V= L/beta*(30*x.^4-60*x.^3+30*x.^2);
            A= L/beta^2*(120*x.^3-180*x.^2+60*x);

    elseif motionlaw==6 % 7th degree poly
            L=L1-L0;
            S = L0+L*(-20*x.^7+70*x.^6-84*x.^5+35*x.^4);
            V= L/beta*(-140*x.^6+420*x.^5-420*x.^4 + 140*x.^3);
            A= L/beta^2*(-840*x.^5+2100*x.^4-1680*x.^3 + 420*x.^2);
    end

    rho = ((R0+S).^2 + V.^2).^(3/2) ./ ((R0+S).^2 + 2*V.^2 - (R0+S).*A);
    rho_min(i) = min(abs(rho));

end

rho_min_beta = rho_min(beta_test*10);

figure;
plot(beta_vec, rho_min);
hold on
plot(beta_test, rho_min(beta_test*10),'ro');
hold off
grid
xlabel('beta (degrees)')
ylabel('rho_{min} (mm)')
end
