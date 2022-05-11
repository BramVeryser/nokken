clear; close all; clc
%% 0-45  (A)

beta_A = 45; %[degree]
start_lift_A = 0; %start lift in mm
end_lift_A = 0; %end lift in mm
motionlaw_A = 1; %dwell

%% 45-120  (B)

beta_B = 75; %[degree]
start_lift_B = 0; %start lift in mm
end_lift_B = 15; %end lift in mm
motionlaw_B = 6;%7th degree poly

%% 120-180  (C)

beta_C = 60; %[degree]
start_lift_C = 15; %start lift in mm
end_lift_C = 30; %end lift in mm
motionlaw_C = 5; %5th degree poly

%% 180-200  (D)

beta_D = 20; %[degree]
start_lift_D = 30; %start lift in mm
end_lift_D = 30; %end lift in mm
motionlaw_D = 1; %dwell


%% 200-280  (E)

beta_E = 80; %[degree]
start_lift_E = 30; %start lift in mm
end_lift_E = 0; %end lift in mm
motionlaw_E = 5; %5th degree poly

%% 280-360   (F)

beta_F = 80; %[degree]
start_lift_F = 0; %start lift in mm
end_lift_F = 0; %end lift in mm
motionlaw_F = 1; %dwell



R0_A = gen_fig_Kloomok_Muffley_1(beta_A,start_lift_A,end_lift_A,motionlaw_A);
R0_B = gen_fig_Kloomok_Muffley_1(beta_B,start_lift_B,end_lift_B,motionlaw_B);
R0_C = gen_fig_Kloomok_Muffley_1(beta_C,start_lift_C,end_lift_C,motionlaw_C);
R0_D = gen_fig_Kloomok_Muffley_1(beta_D,start_lift_D,end_lift_D,motionlaw_D);
R0_E = gen_fig_Kloomok_Muffley_1(beta_E,start_lift_E,end_lift_E,motionlaw_E);
R0_F = gen_fig_Kloomok_Muffley_1(beta_F,start_lift_F,end_lift_F,motionlaw_F);

R0_exact = max([R0_A,R0_B,R0_C,R0_D,R0_E,R0_F])

R0_afgerond = R0_exact + 3.7  %afgerond naar 60

rho_min_A = gen_fig_Kloomok_Muffley_2(R0_afgerond,	,start_lift_A,end_lift_A,motionlaw_A);
rho_min_B = gen_fig_Kloomok_Muffley_2(R0_afgerond,beta_B,start_lift_B,end_lift_B,motionlaw_B);
rho_min_C = gen_fig_Kloomok_Muffley_2(R0_afgerond,beta_C,start_lift_C,end_lift_C,motionlaw_C);
rho_min_D = gen_fig_Kloomok_Muffley_2(R0_afgerond,beta_D,start_lift_D,end_lift_D,motionlaw_D);
rho_min_E = gen_fig_Kloomok_Muffley_2(R0_afgerond,beta_E,start_lift_E,end_lift_E,motionlaw_E);
rho_min_F = gen_fig_Kloomok_Muffley_2(R0_afgerond,beta_F,start_lift_F,end_lift_F,motionlaw_F);

rho_min = min([rho_min_A,rho_min_B,rho_min_C,rho_min_D,rho_min_E,rho_min_F])