clc
clearvars
close all
%% sa = 1 to save figure
sa = 0;
load([pwd '/Eff_density_var_theta.mat'])
%% Parameters

num = 1;

omega  = 2*pi*freq;

theta_i = theta_p(num);

varphi   = a/h;
rho_arit = varphi*rho_1 + (1 - varphi)*rho_2;
rho_geom = (varphi/rho_1 + (1 - varphi)/rho_2)^(-1);

rho_0x = rho_x_eff(num,:);
rho_0z = rho_arit;

%% Plot of effective density rho_x

rho_x_min = 4600;
rho_x_max = 5400;

fts=20;
figure1=figure;
plot1=semilogx(mu_0/lambda_0, rho_x_eff(3,:),...
               mu_0/lambda_0, rho_x_eff(4,:),...
               mu_0/lambda_0, rho_x_eff(5,:),...
               mu_0/lambda_0, rho_arit*ones(length(mu_0),1),...
               mu_0/lambda_0, rho_geom*ones(length(mu_0),1));
xlim([10^(-13) 10^(-4.999999999)])
ylim([rho_x_min rho_x_max])

xlabel('$\mu/\lambda$','fontsize',25,'interpreter','latex')
ylabel('Effective density $\rho_{\textnormal{eff} \parallel}$','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle',':','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[0.2 0.2 0.2]);
set(plot1(5),'Marker','none','Markersize',6,'LineStyle',':','LineWidth',2.0,'Color',[0.2 0.2 0.2]);

legend({'$\rho_{\textnormal{eff} \parallel}$',...
        '$\rho_{a}$','$\rho_{g}$'},...
        'interpreter','latex',...
        'location','east',...
        'fontsize',fts)

set(gca,'fontsize',20,'TickLabelInterpreter','latex','XTick',[1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5])

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    savefig(figure1,'Fig_elas_eff_density_layered_media_var_theta_i_conf.fig')
    print(figure1,'-vector','-dsvg',['Fig_elas_eff_density_layered_media_var_theta_i_conf','.svg'])
end

%% Plotting diference between effective description for different angles
fts=20;
figure1=figure;
plot1=semilogx(mu_0/lambda_0, abs((rho_x_eff(5,:) - rho_x_eff(1,:))/max(rho_x_eff(5,:))), ...
               mu_0/lambda_0, abs((rho_x_eff(5,:) - rho_x_eff(2,:))/max(rho_x_eff(5,:))), ...
               mu_0/lambda_0, abs((rho_x_eff(5,:) - rho_x_eff(3,:))/max(rho_x_eff(5,:))), ...
               mu_0/lambda_0, abs((rho_x_eff(5,:) - rho_x_eff(4,:))/max(rho_x_eff(5,:))));
xlim([10^(-13) 10^(-4.999999999)])

xlabel('$\mu_{0}/\lambda_0$','fontsize',25,'interpreter','latex')
ylabel('Relative error $\eta(\theta_i)$','fontsize',25,'interpreter','latex')

set(plot1(1),'Marker','none','Markersize',6,'LineStyle','-','LineWidth',2.0,'Color',[0 204/255 0]);
set(plot1(2),'Marker','none','Markersize',6,'LineStyle','-.','LineWidth',2.0,'Color',[0 0 204/255]);
set(plot1(3),'Marker','none','Markersize',6,'LineStyle','--','LineWidth',2.0,'Color',[204/255 0 0]);
set(plot1(4),'Marker','none','Markersize',6,'LineStyle',':','LineWidth',2.0,'Color',[204/255 0 204/255]);

legend({'$\eta (\theta_i = \pi/16)$','$\eta (\theta_i = \pi/8)$','$\eta (\theta_i = \pi/6)$','$\eta (\theta_i = \pi/4)$'},...
        'interpreter','latex',...
        'location','northeast',...
        'fontsize',fts)

set(gca,'fontsize',20,'TickLabelInterpreter','latex','XTick',[1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5])

grid on
ax = gca;
ax.XMinorGrid = 'off';

x0=100;
y0=100;
width=750;
height=750;
set(gcf,'position',[x0,y0,width,height])

if sa == 1
    savefig(figure1,'Fig_elas_eff_density_layered_media_var_theta_i_error.fig')
    print(figure1,'-vector','-dsvg',['Fig_elas_eff_density_layered_media_var_theta_i_error','.svg'])
end
