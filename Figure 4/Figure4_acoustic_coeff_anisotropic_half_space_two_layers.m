clc
clearvars
close all
%% sa = 1 to save figures
sa = 0;
%% Physical parameters of layers

freq   = logspace(-4,3,1e4);
omega  = 2*pi*freq;

theta_i = pi/4;

rho_1 = 2500;
K_1   = 2e9;
c_1   = sqrt(K_1/rho_1);

rho_2 = 6000;
K_2   = 1.75*K_1;
c_2   = sqrt(K_2/rho_2);

h = 0.1;
a = h/4;
b = h - a;

varphi = a/h;

rho_0z = 1.2;
rho_0x_g = 10*rho_0z;
rho_0x_l = 0.1*rho_0z;

K_0   = 1.4e5;
c_0   = sqrt(K_0/rho_0z);
k_0   = omega/c_0;

%% Reflection and transmission for 2 layers

[R_g, T_g, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Acoustic_reflection_transmission_anisotropic_2_layers(freq, rho_0x_g, rho_0z, rho_1, rho_2, K_0, K_1, K_2, h, varphi, theta_i); 

[R_l, T_l, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Acoustic_reflection_transmission_anisotropic_2_layers(freq, rho_0x_l, rho_0z, rho_1, rho_2, K_0, K_1, K_2, h, varphi, theta_i); 

[R_iso, T_iso, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = Acoustic_reflection_transmission_2_layers(freq, rho_0z, rho_1, rho_2, K_0, K_1, K_2, h, varphi, theta_i); 

%% Plot of reflection coefficient (real and imaginary parts)

fts=20;
figure1=figure;
plot1=semilogx(k_0*h, real(R_l),...
               k_0*h, -imag(R_l),...
               k_0*h, real(R_iso),...
               k_0*h, -imag(R_iso),...
               k_0*h, real(R_g),...
               k_0*h, -imag(R_g));
xlim([1e-5 1e-1])
ylim([0 1])

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

set(gca,'fontsize',18,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

xlabel('$k_0 (\theta_i) h$','fontsize',25,'interpreter','latex')
ylabel('Reflection coefficient','fontsize',25,'interpreter','latex')

legend({'$\mathrm{Re}\{R\}$ -- $\mathcal{R}_\rho = 0.1$','$-\mathrm{Im}\{R\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$\mathrm{Re}\{R\}$ -- isotropic','$-\mathrm{Im}\{R\}$ -- isotropic',...
        '$\mathrm{Re}\{R\}$ -- $\mathcal{R}_\rho = 10$','$-\mathrm{Im}\{R\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northwest',...
        'fontsize',fts)

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_R_acous_aniso_layered_media','.svg'])
end

%% Plot of transmission coefficient (real and imaginary parts)
fts=20;
figure1=figure;
plot1=semilogx(k_0*h, real(T_l),...
               k_0*h, imag(T_l),...
               k_0*h, real(T_iso),...
               k_0*h, imag(T_iso),...
               k_0*h, real(T_g),...
               k_0*h, imag(T_g));
xlim([1e-5 1e-1])
ylim([0 1])

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 170/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[77/255 225/255 1]);
set(plot1(6),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[128/255 85/255 0]);

set(gca,'fontsize',18,'xtick',[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1],'TickLabelInterpreter','latex')

xlabel('$k_0 (\theta_i) h$','fontsize',25,'interpreter','latex')
ylabel('Transmission coefficient','fontsize',25,'interpreter','latex')

legend({'$\mathrm{Re}\{T\}$ -- $\mathcal{R}_\rho = 0.1$','$\mathrm{Im}\{T\}$ -- $\mathcal{R}_\rho = 0.1$',...
        '$\mathrm{Re}\{T\}$ -- isotropic','$\mathrm{Im}\{T\}$ -- isotropic',...
        '$\mathrm{Re}\{T\}$ -- $\mathcal{R}_\rho = 10$','$\mathrm{Im}\{T\}$ -- $\mathcal{R}_\rho = 10$'},...
        'interpreter','latex',...
        'location','northeast',...
        'fontsize',fts)

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    print(figure1,'-vector','-dsvg',['Fig_T_acous_aniso_layered_media','.svg'])
end
