clc
clearvars
close all
%% sa = 1 to save figures
sa = 0;
%% Physical parameters of material

freq   = logspace(-2,3,1e2);
omega  = 2*pi*freq;

rho_0z = 1.2;
rho_0x_g = 10*rho_0z;
rho_0x_l = 0.1*rho_0z;
rho_1 = 2500;
rho_2 = 6000;

mu_0     = 0.75e8;
lambda_0 = 1e9;

mu_1     = 5e9;
lambda_1 = 20e9;

mu_2     = 0.5e9;
lambda_2 = 5e9;

theta_p = pi/4;

a = 0.025;
h = 0.1;

for uf = 1:length(omega)
    
    [R_p_g(uf), R_s_g(uf), T_p_g(uf), T_s_g(uf), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        Reflection_transmission_elasticity_aniso_layered_media(omega(uf), mu_0, lambda_0, rho_0x_g, rho_0z, mu_1, lambda_1, rho_1, mu_2, lambda_2, rho_2, a, h, theta_p);

    [R_p_iso(uf), R_s_iso(uf), T_p_iso(uf), T_s_iso(uf), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        Reflection_transmission_elasticity_layered_media(omega(uf), mu_0, lambda_0, rho_0z, mu_1, lambda_1, rho_1, mu_2, lambda_2, rho_2, a, h, theta_p);

    [R_p_l(uf), R_s_l(uf), T_p_l(uf), T_s_l(uf), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        Reflection_transmission_elasticity_aniso_layered_media(omega(uf), mu_0, lambda_0, rho_0x_l, rho_0z, mu_1, lambda_1, rho_1, mu_2, lambda_2, rho_2, a, h, theta_p);
end

c_p0 = sqrt((lambda_0 + 2*mu_0)/rho_0z);
k_p0 = omega/c_p0;

%% Plot of reflection coefficient for p-waves
fts=20;
figure1=figure;
plot1=semilogx(k_p0*h, real(R_p_l),...
               k_p0*h, -imag(R_p_l),...
               k_p0*h, real(R_p_iso),...
               k_p0*h, -imag(R_p_iso),...
               k_p0*h, real(R_p_g),...
               k_p0*h, -imag(R_p_g));
xlim([1e-6 1e-2])
ylim([-0.8 1])

xlabel('$k_{0p} h$','fontsize',25,'interpreter','latex')
ylabel('Reflection coefficients','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

legend({'$\mathrm{Re}\{R_p\}$ -- $\mathcal{R}_\rho = 0.1$','$-\mathrm{Im}\{R_p\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$\mathrm{Re}\{R_p\}$ -- isotropic','$-\mathrm{Im}\{R_p\}$ -- isotropic',...
        '$\mathrm{Re}\{R_p\}$ -- $\mathcal{R}_\rho = 10$','$-\mathrm{Im}\{R_p\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northwest',...
        'fontsize',fts)

set(gca,'fontsize',20,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_R_p_elas_aniso_layered_media','.svg'])
end

%% Plot of reflection coefficient for s-waves
fts=20;
figure1=figure;
plot1=semilogx(k_p0*h, -real(R_s_l),...
               k_p0*h, imag(R_s_l),...
               k_p0*h, -real(R_s_iso),...
               k_p0*h, imag(R_s_iso),...
               k_p0*h, -real(R_s_g),...
               k_p0*h, imag(R_s_g));
xlim([1e-6 1e-2])
ylim([0 1.4])

xlabel('$k_{0p} h$','fontsize',25,'interpreter','latex')
ylabel('Reflection coefficients','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

legend({'$-\mathrm{Re}\{R_s\}$ -- $\mathcal{R}_\rho = 0.1$','$\mathrm{Im}\{R_s\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$-\mathrm{Re}\{R_s\}$ -- isotropic','$\mathrm{Im}\{R_s\}$ -- isotropic',...
        '$-\mathrm{Re}\{R_s\}$ -- $\mathcal{R}_\rho = 10$','$\mathrm{Im}\{R_s\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northwest',...
        'fontsize',fts)

set(gca,'fontsize',20,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_R_s_elas_aniso_layered_media','.svg'])
end

%% Plot of transmission coefficient for p-waves
fts=20;
figure1=figure;
plot1=semilogx(k_p0*h, real(T_p_l),...
               k_p0*h, imag(T_p_l),...
               k_p0*h, real(T_p_iso),...
               k_p0*h, imag(T_p_iso),...
               k_p0*h, real(T_p_g),...
               k_p0*h, imag(T_p_g));
xlim([1e-6 1e-2])
ylim([0 1])

xlabel('$k_{0p} h$','fontsize',25,'interpreter','latex')
ylabel('Transmission coefficients','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

legend({'$\mathrm{Re}\{T_p\}$ -- $\mathcal{R}_\rho = 0.1$','$\mathrm{Im}\{T_p\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$\mathrm{Re}\{T_p\}$ -- isotropic','$\mathrm{Im}\{T_p\}$ -- isotropic',...
        '$\mathrm{Re}\{T_p\}$ -- $\mathcal{R}_\rho = 10$','$\mathrm{Im}\{T_p\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northeast',...
        'fontsize',fts)

set(gca,'fontsize',20,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_T_p_elas_aniso_layered_media','.svg'])
end

%% Plot of transmission coefficient for s-waves
fts=20;
figure1=figure;
plot1=semilogx(k_p0*h, real(T_s_l),...
               k_p0*h, imag(T_s_l),...
               k_p0*h, real(T_s_iso),...
               k_p0*h, imag(T_s_iso),...
               k_p0*h, real(T_s_g),...
               k_p0*h, imag(T_s_g));
xlim([1e-6 1e-2])
ylim([-0.4 0.6])

xlabel('$k_{0p} h$','fontsize',25,'interpreter','latex')
ylabel('Transmission coefficients','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

legend({'$\mathrm{Re}\{T_s\}$ -- $\mathcal{R}_\rho = 0.1$','$\mathrm{Im}\{T_s\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$\mathrm{Re}\{T_s\}$ -- isotropic','$\mathrm{Im}\{T_s\}$ -- isotropic',...
        '$\mathrm{Re}\{T_s\}$ -- $\mathcal{R}_\rho = 10$','$\mathrm{Im}\{T_s\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northwest',...
        'fontsize',fts)

set(gca,'fontsize',20,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_T_s_elas_aniso_layered_media','.svg'])
end

%% Acoustic limit

omega = 100;
freq = omega/(2*pi);

varphi = a/h;

mu_exp_min = 0;
mu_exp_max = 9;
mu_1       = logspace(mu_exp_min,mu_exp_max,1e2);
mu_2       = mu_1;
mu_0       = mu_1;

m = 0;

for um = 1:length(mu_1)
    [R_p(um), ~, T_p(um), ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        Reflection_transmission_elasticity_aniso_layered_media(omega, mu_0(um), lambda_0, rho_0x_l, rho_0z, mu_1(um), lambda_1, rho_1, mu_2(um), lambda_2, rho_2, a, h, theta_p);
end

[R, T, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Acoustic_reflection_transmission_anisotropic_2_layers(freq, rho_0x_l, rho_0z, rho_1, rho_2, lambda_0, lambda_1, lambda_2, h, varphi, theta_p); 

%% Plot of acoustic limit for p-waves
fts=20;
figure1=figure;
plot1=semilogx(mu_0/lambda_0, abs(R_p),...
               mu_0/lambda_0, abs(T_p),...
               mu_0/lambda_0, abs(R)*ones(length(mu_1),1),...
               mu_0/lambda_0, abs(T)*ones(length(mu_1),1));
xlim([1e-8 1])

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle',':','LineWidth',2.0,'Color',[0 0 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[0 0 0]);

xlabel('$\mu/\lambda_0$','fontsize',fts,'interpreter','latex')
ylabel('Reflection and transmission coefficients','fontsize',fts,'interpreter','latex')

legend({'$|R_{p}|$','$|T_{p}|$', '$|R|$ -- acoustics', '$|T|$ -- acoustics'},...
        'interpreter','latex',...
        'location','west',...
        'fontsize',fts)

set(gca,'fontsize',20,'TickLabelInterpreter','latex')

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    saveas(figure1,'Fig_R_T_elas_acoustic_limit_aniso_layered_media','svg')
end