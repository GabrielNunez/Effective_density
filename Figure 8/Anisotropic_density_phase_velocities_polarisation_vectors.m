function [c, A_p, A_s] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda, mu, rho_x, rho_z, theta)

% if mu < 1e0
%     mu = vpa(mu);
% end

Gamma_11 = mu*cos(theta)^2 + (lambda + 2*mu)*sin(theta)^2;
Gamma_33 = (lambda + 2*mu)*cos(theta)^2 + mu*sin(theta)^2;
Gamma_13 = (lambda + mu)*cos(theta)*sin(theta);

R_rho = rho_x/rho_z;

c_p = sqrt((Gamma_33*R_rho + Gamma_11 + sqrt((Gamma_33*R_rho - Gamma_11).^2 + 4*Gamma_13.^2*R_rho))./(2*rho_x));
c_s = sqrt((Gamma_33*R_rho + Gamma_11 - sqrt((Gamma_33*R_rho - Gamma_11).^2 + 4*Gamma_13.^2*R_rho))./(2*rho_x));

% c_p = double(c_p);
% c_s = double(c_s);

if theta == 0
    A_p1 = 0;
    A_p3 = 1;
    A_s1 = 1;
    A_s3 = 0;
elseif theta == pi/2
    A_p1 = 1;
    A_p3 = 0;
    A_s1 = 0;
    A_s3 = 1;
else
    Sigma = (Gamma_33*rho_x - Gamma_11*rho_z)./(2*Gamma_13*rho_z);
    A_p1 =  1./sqrt(1 + (Sigma + sqrt(Sigma.^2 + rho_x./rho_z)).^2);
    A_p3 =  (Sigma + sqrt(Sigma.^2 + rho_x./rho_z))./sqrt(1 + (Sigma + sqrt(Sigma.^2 + rho_x./rho_z)).^2);
    A_s1 =  - 1./sqrt(1 + (Sigma - sqrt(Sigma.^2 + rho_x./rho_z)).^2);
    A_s3 =  (- Sigma + sqrt(Sigma.^2 + rho_x./rho_z))./sqrt(1 + (Sigma - sqrt(Sigma.^2 + rho_x./rho_z)).^2);
end

% A_p1 = double(A_p1);
% A_p3 = double(A_p3);
% A_s1 = double(A_s1);
% A_s3 = double(A_s3);

c     = [c_p; c_s];
A_p   = [A_p1; A_p3];
A_s   = [A_s1; A_s3];