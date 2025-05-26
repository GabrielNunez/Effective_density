function [R_p, R_s, T_p, T_s, R_p1, R_s1, T_p1, T_s1, R_p2, R_s2, T_p2, T_s2, J_pr, J_sr, J_pt, J_st] = ...
    Reflection_transmission_elasticity_aniso_layered_media_direct(omega, mu_0, lambda_0, rho_0x, rho_0z, mu_1, lambda_1, rho_1, mu_2, lambda_2, rho_2, a, h, theta_i)

[c_1, P_p, ~] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda_0, mu_0, rho_0x, rho_0z, theta_i);

c_p0  = c_1(1);
P_p0x = P_p(1);
P_p0z = P_p(2);

options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-16);

if theta_i == 0
    theta_s = 0;
else
    if mu_0/lambda_0 > 1e-8
        x = fsolve(@(theta_s)(sqrt((((lambda_0 + 2*mu_0)*cos(theta_s).^2 + mu_0*sin(theta_s).^2)*rho_0x + (mu_0*cos(theta_s).^2 + (lambda_0 + 2*mu_0)*sin(theta_s).^2)*rho_0z - sqrt((((lambda_0 + 2*mu_0)*cos(theta_s).^2 + mu_0*sin(theta_s).^2)*rho_0x - (mu_0*cos(theta_s).^2 + (lambda_0 + 2*mu_0)*sin(theta_s).^2)*rho_0z).^2 + 4*((lambda_0 + mu_0)*cos(theta_s).*sin(theta_s)).^2*rho_0x*rho_0z))./(2*rho_0x*rho_0z))/sin(theta_s) - c_p0./sin(theta_i)),1e-2,options);
        theta_s = x;
    else
        theta_s = 0;
    end
end

[c_1, ~, P_s] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda_0, mu_0, rho_0x, rho_0z, theta_s);
c_s0  = c_1(2);
P_s0x = P_s(1);
P_s0z = P_s(2);

c_p1 = sqrt((lambda_1 + 2*mu_1)/rho_1);
c_s1 = sqrt(mu_1/rho_1);

c_p2 = sqrt((lambda_2 + 2*mu_2)/rho_2);
c_s2 = sqrt(mu_2/rho_2);

theta_p1 = asin(c_p1/c_p0*sin(theta_i));
theta_s1 = asin(c_s1/c_p0*sin(theta_i));
theta_p2 = asin(c_p2/c_p0*sin(theta_i));
theta_s2 = asin(c_s2/c_p0*sin(theta_i));

P_p1x    = sin(theta_p1);
P_p1z    = cos(theta_p1);
P_s1x    = - cos(theta_s1);
P_s1z    = sin(theta_s1);
P_p2x    = sin(theta_p2);
P_p2z    = cos(theta_p2);
P_s2x    = - cos(theta_s2);
P_s2z    = sin(theta_s2);

k_p0 = omega/c_p0;
k_s0 = omega/c_s0;
k_p1 = omega/c_p1;
k_s1 = omega/c_s1;
k_p2 = omega/c_p2;
k_s2 = omega/c_s2;

A = zeros(12,12);
b = zeros(12,1);

A(1,1) =   P_p0x;
A(1,2) = - P_s0x;
A(1,3) = - P_p1x;
A(1,4) = - P_p1x;
A(1,5) = - P_s1x;
A(1,6) =   P_s1x;

A(2,1) = - P_p0z;
A(2,2) =   P_s0z;
A(2,3) = - P_p1z;
A(2,4) =   P_p1z;
A(2,5) = - P_s1z;
A(2,6) = - P_s1z;

A(3,1) =   k_p0*(P_p0x*lambda_0*sin(theta_i) + P_p0z*(lambda_0 + 2*mu_0)*cos(theta_i));
A(3,2) = - k_s0*(P_s0x*lambda_0*sin(theta_s) + P_s0z*(lambda_0 + 2*mu_0)*cos(theta_s));
A(3,3) = - k_p1*(P_p1x*lambda_1*sin(theta_p1) + P_p1z*(lambda_1 + 2*mu_1)*cos(theta_p1));
A(3,4) = - k_p1*(P_p1x*lambda_1*sin(theta_p1) + P_p1z*(lambda_1 + 2*mu_1)*cos(theta_p1));
A(3,5) = - k_s1*(P_s1x*lambda_1*sin(theta_s1) + P_s1z*(lambda_1 + 2*mu_1)*cos(theta_s1));
A(3,6) =   k_s1*(P_s1x*lambda_1*sin(theta_s1) + P_s1z*(lambda_1 + 2*mu_1)*cos(theta_s1));

A(4,1) = - mu_0*k_p0*(P_p0x*cos(theta_i) + P_p0z*sin(theta_i));
A(4,2) = mu_0*k_s0*(P_s0x*cos(theta_s) + P_s0z*sin(theta_s));
A(4,3) = - mu_1*k_p1*(P_p1x*cos(theta_p1) + P_p1z*sin(theta_p1));
A(4,4) = mu_1*k_p1*(P_p1x*cos(theta_p1) + P_p1z*sin(theta_p1));
A(4,5) = - mu_1*k_s1*(P_s1x*cos(theta_s1) + P_s1z*sin(theta_s1));
A(4,6) = - mu_1*k_s1*(P_s1x*cos(theta_s1) + P_s1z*sin(theta_s1));

A(5,3)  = P_p1x*exp(1i*k_p1*a*cos(theta_p1));
A(5,4)  = P_p1x*exp(- 1i*k_p1*a*cos(theta_p1));
A(5,5)  = P_s1x*exp(1i*k_s1*a*cos(theta_s1));
A(5,6)  = - P_s1x*exp(- 1i*k_s1*a*cos(theta_s1));
A(5,7)  = - P_p2x*exp(1i*k_p2*a*cos(theta_p2));
A(5,8)  = - P_p2x*exp(- 1i*k_p2*a*cos(theta_p2));
A(5,9)  = - P_s2x*exp(1i*k_s2*a*cos(theta_s2));
A(5,10) = P_s2x*exp(- 1i*k_s2*a*cos(theta_s2));

A(6,3)  = P_p1z*exp(1i*k_p1*a*cos(theta_p1));
A(6,4)  = - P_p1z*exp(- 1i*k_p1*a*cos(theta_p1));
A(6,5)  = P_s1z*exp(1i*k_s1*a*cos(theta_s1));
A(6,6)  = P_s1z*exp(- 1i*k_s1*a*cos(theta_s1));
A(6,7)  = - P_p2z*exp(1i*k_p2*a*cos(theta_p2));
A(6,8)  = P_p2z*exp(- 1i*k_p2*a*cos(theta_p2));
A(6,9)  = - P_s2z*exp(1i*k_s2*a*cos(theta_s2));
A(6,10) = - P_s2z*exp(- 1i*k_s2*a*cos(theta_s2));

A(7,3)  = k_p1*(P_p1x*lambda_1*sin(theta_p1) + P_p1z*(lambda_1 + 2*mu_1)*cos(theta_p1))*exp(1i*k_p1*a*cos(theta_p1));
A(7,4)  = k_p1*(P_p1x*lambda_1*sin(theta_p1) + P_p1z*(lambda_1 + 2*mu_1)*cos(theta_p1))*exp(- 1i*k_p1*a*cos(theta_p1));
A(7,5)  = k_s1*(P_s1x*lambda_1*sin(theta_s1) + P_s1z*(lambda_1 + 2*mu_1)*cos(theta_s1))*exp(1i*k_s1*a*cos(theta_s1));
A(7,6)  = - k_s1*(P_s1x*lambda_1*sin(theta_s1) + P_s1z*(lambda_1 + 2*mu_1)*cos(theta_s1))*exp(- 1i*k_s1*a*cos(theta_s1));
A(7,7)  = - k_p2*(P_p2x*lambda_2*sin(theta_p2) + P_p2z*(lambda_2 + 2*mu_2)*cos(theta_p2))*exp(1i*k_p2*a*cos(theta_p2));
A(7,8)  = - k_p2*(P_p2x*lambda_2*sin(theta_p2) + P_p2z*(lambda_2 + 2*mu_2)*cos(theta_p2))*exp(- 1i*kp2*a*cos(theta_p2));
A(7,9)  = - k_s2*(P_s2x*lambda_2*sin(theta_s2) + P_s2z*(lambda_2 + 2*mu_2)*cos(theta_s2))*exp(1i*k_s2*a*cos(theta_s2));
A(7,10) = k_s2*(P_s2x*lambda_2*sin(theta_s2) + P_s2z*(lambda_2 + 2*mu_2)*cos(theta_s2))*exp(- 1i*k_s2*a*cos(theta_s2));

A(8,3)  = mu_1*k_p1*(P_p1x*cos(theta_p1) + P_p1z*sin(theta_p1))*exp(1i*k_p1*a*cos(theta_p1));
A(8,4)  = - mu_1*k_p1*(P_p1x*cos(theta_p1) + P_p1z*sin(theta_p1))*exp(- 1i*k_p1*a*cos(theta_p1));
A(8,5)  = mu_1*k_s1*(P_s1x*cos(theta_s1) + P_s1z*sin(theta_s1))*exp(1i*k_s1*a*cos(theta_s1));
A(8,6)  = mu_1*k_s1*(P_s1x*cos(theta_s1) + P_s1z*sin(theta_s1))*exp(- 1i*k_s1*a*cos(theta_s1));
A(8,7)  = - mu_2*k_p2*(P_p2x*cos(theta_p2) + P_p2z*sin(theta_p2))*exp(1i*k_p2*a*cos(theta_p2));
A(8,8)  = mu_2*k_p2*(P_p2x*cos(theta_p2) + P_p2z*sin(theta_p2))*exp(- 1i*k_p2*a*cos(theta_p2));
A(8,9)  = - mu_2*k_s2*(P_s2x*cos(theta_s2) + P_s2z*sin(theta_s2))*exp(1i*k_s2*a*cos(theta_s2));
A(8,10) = - mu_2*k_s2*(P_s2x*cos(theta_s2) + P_s2z*sin(theta_s2))*exp(- 1i*k_s2*a*cos(theta_s2));

A(9,7)  = P_p2x*exp(1i*k_p2*h*cos(theta_p2));
A(9,8)  = P_p2x*exp(- 1i*k_p2*h*cos(theta_p2));
A(9,9)  = P_s2x*exp(1i*k_s2*h*cos(theta_s2));
A(9,10) = - P_s2x*exp(- 1i*k_s2*h*cos(theta_s2));
A(9,11) = - P_p0x*exp(1i*k_p0*h*cos(theta_i));
A(9,12) = - P_s0x*exp(1i*k_s0*h*cos(theta_s));

A(10,7)  = P_p2z*exp(1i*k_p2*h*cos(theta_p2));
A(10,8)  = - P_p2z*exp(- 1i*k_p2*h*cos(theta_p2));
A(10,9)  = P_s2z*exp(1i*k_s2*h*cos(theta_s2));
A(10,10) = P_s2z*exp(- 1i*k_s2*h*cos(theta_s2));
A(10,11) = - P_p0z*exp(1i*k_p0*h*cos(theta_i));
A(10,12) = - P_s0z*exp(1i*k_s0*h*cos(theta_s));

A(11,7)  = k_p2*(P_p2x*lambda_2*sin(theta_p2) + P_p2z*(lambda_2 + 2*mu_2)*cos(theta_p2))*exp(1i*k_p2*h*cos(theta_p2));
A(11,8)  = k_p2*(P_p2x*lambda_2*sin(theta_p2) + P_p2z*(lambda_2 + 2*mu_2)*cos(theta_p2))*exp(- 1i*k_p2*h*cos(theta_p2));
A(11,9)  = k_s2*(P_s2x*lambda_2*sin(theta_s2) + P_s2z*(lambda_2 + 2*mu_2)*cos(theta_s2))*exp(1i*k_s2*h*cos(theta_s2));
A(11,10) = - k_s2*(P_s2x*lambda_2*sin(theta_s2) + P_s2z*(lambda_2 + 2*mu_2)*cos(theta_s2))*exp(- 1i*k_s2*h*cos(theta_s2));
A(11,11) = - k_p0*(P_p0x*lambda_0*sin(theta_i) + P_p0z*(lambda_0 + 2*mu_0)*cos(theta_i))*exp(1i*k_p0*h*cos(theta_i));
A(11,12) = - k_s0*(P_s0x*lambda_0*sin(theta_s) + P_s0z*(lambda_0 + 2*mu_0)*cos(theta_s))*exp(1i*k_s0*h*cos(theta_s));

A(12,7)  = mu_2*k_p2*(P_p2x*cos(theta_p2) + P_p2z*sin(theta_p2))*exp(1i*k_p2*h*cos(theta_p2));
A(12,8)  = - mu_2*k_p2*(P_p2x*cos(theta_p2) + P_p2z*sin(theta_p2))*exp(- 1i*k_p2*h*cos(theta_p2));
A(12,9)  = mu_2*k_s2*(P_s2x*cos(theta_s2) + P_s2z*sin(theta_s2))*exp(1i*k_s2*h*cos(theta_s2));
A(12,10) = mu_2*k_s2*(P_s2x*cos(theta_s2) + P_s2z*sin(theta_s2))*exp(- 1i*k_s2*h*cos(theta_s2));
A(12,11) = - mu_0*k_p0*(P_p0x*cos(theta_i) + P_p0z*sin(theta_i))*exp(1i*k_p0*h*cos(theta_i));
A(12,12) = - mu_0*k_s0*(P_s0x*cos(theta_s) + P_s0z*sin(theta_s))*exp(1i*k_s0*h*cos(theta_s));

b(1,1) = - P_p0x;
b(2,1) = - P_p0z;
b(3,1) = - k_p0*(P_p0x*lambda_0*sin(theta_i) + P_p0z*(lambda_0 + 2*mu_0)*cos(theta_i));
b(4,1) = - mu_0*k_p0*(P_p0x*cos(theta_i) + P_p0z*sin(theta_i));

X(:)    = linsolve(A,b);
R_p   = X(1);
R_s   = X(2);
T_p1  = X(3);
R_p1  = X(4);
T_s1  = X(5);
R_s1  = X(6);
T_p2  = X(7);
R_p2  = X(8);
T_s2  = X(9);
R_s2  = X(10);
T_p   = X(11);
T_s   = X(12);

J_pr = abs(R_p).^2;
J_sr = abs(R_s).^2*c_p0/c_s0*(mu_0*cos(theta_s)*P_s0x^2 + (lambda_0 + mu_0)*sin(theta_s)*P_s0x*P_s0z + (lambda_0 + 2*mu_0)*cos(theta_s)*P_s0z^2)/(mu_0*cos(theta_i)*P_p0x^2 + (lambda_0 + mu_0)*sin(theta_i)*P_p0x*P_p0z + (lambda_0 + 2*mu_0)*cos(theta_i)*P_p0z^2);
J_pt = abs(T_p).^2;
J_st = abs(T_s).^2*c_p0/c_s0*(mu_0*cos(theta_s)*P_s0x^2 + (lambda_0 + mu_0)*sin(theta_s)*P_s0x*P_s0z + (lambda_0 + 2*mu_0)*cos(theta_s)*P_s0z^2)/(mu_0*cos(theta_i)*P_p0x^2 + (lambda_0 + mu_0)*sin(theta_i)*P_p0x*P_p0z + (lambda_0 + 2*mu_0)*cos(theta_i)*P_p0z^2);
