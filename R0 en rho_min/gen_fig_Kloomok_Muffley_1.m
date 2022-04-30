%this script is a replacement for the diagram on slide 33 of chapter 8
%this script will generate graphs that can be used for all motion laws

function [R_0_30] = gen_fig_Kloomok_Muffley_1(beta,start_lift,end_lift,motionlaw)

%% input
% beta = 75; %[degree]
% start_lift = 0; %start lift in mm
% end_lift = 15; %end lift in mm
% motionlaw = 4; 

%% calculations
R0_vec = 0:0.1:100;
alpha_max = 0*R0_vec;

L0 = start_lift;
L1 = end_lift;
beta = beta*pi/180;
for i = 1:length(R0_vec)
    R0 = R0_vec(i);
    x = 0:0.01:1;
    
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

    alpha = atan2(V,R0+S);
    alpha_max(i) = max(abs(alpha))*180/pi;

end

R_0_30 = find(abs(alpha_max-30) == min(abs(alpha_max-30)),1)/10;

figure;
plot(R0_vec, alpha_max);
grid
xlabel('R_0 (mm)')
ylabel('alpha_{max} (degree)')
end
