clc
clearvars
close all
%% sa = 1 for saving effective density values
sa = 0;
%% Parameters

freq   = 1;
omega  = 2*pi*freq;

rho_1 = 3000;
rho_2 = 6000;

mu_1     = logspace(-4,6,1e3);
lambda_1 = 1e11;

mu_2     = mu_1;
lambda_2 = lambda_1;

mu_0     = mu_1;
lambda_0 = lambda_1;

theta_p = [pi/16,pi/8,pi/6,pi/4,pi/3];

a = 0.025;
h = 0.1;

varphi   = a/h;
rho_arit = varphi*rho_1 + (1 - varphi)*rho_2;
rho_geom = (varphi/rho_1 + (1 - varphi)/rho_2)^(-1);

rho_x_lim = 200;

rho_x_min = 4600;
rho_x_max = 5400;

rho_0z = rho_arit;

%% Loop for determining effective density rho_x

step = 10;
bs   = step;

rho_0x_g = 4800;

n = 0;
m = 0;

max_n = 4;

for ui = 1:length(theta_p)
    for ua = 1:length(mu_1)

        rho_0x_g_pl = rho_0x_g + bs;
        rho_0x_g_mi = rho_0x_g - bs;

        [R_p_g,~, T_p_g,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                    Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
            
        F_g = abs(R_p_g) + abs(1 - T_p_g);

        [R_p_g_pl,~, T_p_g_pl,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_pl, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
        
        F_g_pl = abs(R_p_g_pl) + abs(1 - T_p_g_pl);

        [R_p_g_mi,~, T_p_g_mi,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_mi, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
        
        F_g_mi = abs(R_p_g_mi) + abs(1 - T_p_g_mi);

        F_min = min([F_g,F_g_pl,F_g_mi]);

        while n <= max_n
            
            while F_min ~= F_g

                if F_min == F_g_pl
                    rho_0x_g = rho_0x_g + bs;
                    rho_0x_g_pl = rho_0x_g + bs;
                    rho_0x_g_mi = rho_0x_g - bs;
            
                    F_g_mi = F_g;
                    F_g    = F_g_pl;
            
                    [R_p_g_pl,~, T_p_g_pl,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                            Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_pl, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
                    
                    F_g_pl = abs(R_p_g_pl) + abs(1 - T_p_g_pl);

                    F_min = min([F_g,F_g_pl,F_g_mi]);
                elseif F_min == F_g_mi
                    rho_0x_g = rho_0x_g - bs;
                    rho_0x_g_pl = rho_0x_g + bs;
                    rho_0x_g_mi = rho_0x_g - bs;

                    F_g_pl = F_g;
                    F_g    = F_g_mi;
            
                    [R_p_g_mi,~, T_p_g_mi,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                            Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_mi, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
                    
                    F_g_mi = abs(R_p_g_mi) + abs(1 - T_p_g_mi);
            
                    F_min = min([F_g,F_g_pl,F_g_mi]);
                end
            end

        bs = bs/10;

        rho_0x_g_pl = rho_0x_g + bs;
        rho_0x_g_mi = rho_0x_g - bs;

        [R_p_g_pl,~, T_p_g_pl,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_pl, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
        
        F_g_pl = abs(R_p_g_pl) + abs(1 - T_p_g_pl);

        [R_p_g_mi,~, T_p_g_mi,~,~,~,~,~,~,~,~,~,~,~,~,~] = ...
                Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(ua), lambda_0, rho_0x_g_mi, rho_0z, mu_1(ua), lambda_1, rho_1, mu_2(ua), lambda_2, rho_2, a, h, theta_p(ui));
        
        F_g_mi = abs(R_p_g_mi) + abs(1 - T_p_g_mi);

        F_min = min([F_g,F_g_pl,F_g_mi]);
        
        n = n + 1; 
        end  

        rho_x_eff(ui,ua) = rho_0x_g;

        bs = step;

        n = 0;
        m = m + 1

    end
end

rho_z_eff = rho_arit;

%% Saving .mat file with effective density values
if sa == 1
    save([pwd '/Eff_density_files/Eff_density_var_theta_new_alg.mat'],'freq','mu_0','lambda_0','rho_1','mu_1','lambda_1','rho_2','mu_2','lambda_2','theta_p','a','h','rho_x_eff','rho_z_eff')
end
