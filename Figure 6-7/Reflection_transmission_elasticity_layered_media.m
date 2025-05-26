function [R_p, R_s, T_p, T_s, R_p1, R_s1, T_p1, T_s1, R_p2, R_s2, T_p2, T_s2, J_pr, J_sr, J_pt, J_st] = ...
    Reflection_transmission_elasticity_layered_media_direct(omega, mu_0, lambda_0, rho_0, mu_1, lambda_1, rho_1, mu_2, lambda_2, rho_2, a, h, theta_p)

c_p0 = sqrt((lambda_0 + 2*mu_0)/rho_0);
c_s0 = sqrt(mu_0/rho_0);
c_p1 = sqrt((lambda_1 + 2*mu_1)/rho_1);
c_s1 = sqrt(mu_1/rho_1);
c_p2 = sqrt((lambda_2 + 2*mu_2)/rho_2);
c_s2 = sqrt(mu_2/rho_2);

k_p0 = omega/c_p0;
k_s0 = omega/c_s0;
k_p1 = omega/c_p1;
k_s1 = omega/c_s1;
k_p2 = omega/c_p2;
k_s2 = omega/c_s2;

theta_s   = asin(c_s0/c_p0*sin(theta_p));
theta_p1  = asin(c_p1/c_p0*sin(theta_p));
theta_s1  = asin(c_s1/c_p0*sin(theta_p));
theta_p2  = asin(c_p2/c_p0*sin(theta_p));
theta_s2  = asin(c_s2/c_p0*sin(theta_p));

k_p0x = k_p0*sin(theta_p);
k_s0z = k_s0*cos(theta_s);
k_p1z = k_p1*cos(theta_p1);
k_s1z = k_s1*cos(theta_s1);
k_p2z = k_p2*cos(theta_p2);
k_s2z = k_s2*cos(theta_s2);

r_1 = cot(theta_p);
r_2 = k_s0z/k_p0x;
r_3 = k_p1z/k_p0x;
r_4 = k_s1z/k_p0x;
r_5 = k_p2z/k_p0x;
r_6 = k_s2z/k_p0x;

A = zeros(12,12);
b = zeros(12,1);

A(1,1) = 1;
A(1,2) = r_2;
A(1,3) = - 1;
A(1,4) = - 1;
A(1,5) = r_4;
A(1,6) = - r_4;

A(2,1) = - r_1;
A(2,2) = 1;
A(2,3) = - r_3;
A(2,4) = r_3;
A(2,5) = - 1;
A(2,6) = - 1;

A(3,1) = - (lambda_0 + r_1^2*(lambda_0 + 2*mu_0));
A(3,2) = 2*mu_0*r_2;
A(3,3) = (lambda_1 + r_3^2*(lambda_1 + 2*mu_1));
A(3,4) = (lambda_1 + r_3^2*(lambda_1 + 2*mu_1));
A(3,5) = 2*mu_1*r_4;
A(3,6) = - 2*mu_1*r_4;

A(4,1) = 2*mu_0*r_1;
A(4,2) = mu_0*(r_2^2 - 1);
A(4,3) = 2*mu_1*r_3;
A(4,4) = - 2*mu_1*r_3;
A(4,5) = - mu_1*(r_4^2 - 1);
A(4,6) = - mu_1*(r_4^2 - 1);

A(5,3)  = exp(1i*k_p1*a*cos(theta_p1));
A(5,4)  = exp(-1i*k_p1*a*cos(theta_p1));
A(5,5)  = - r_4*exp(1i*k_s1*a*cos(theta_s1));
A(5,6)  = r_4*exp(-1i*k_s1*a*cos(theta_s1));
A(5,7)  = - exp(1i*k_p2*a*cos(theta_p2));
A(5,8)  = - exp(-1i*k_p2*a*cos(theta_p2));
A(5,9)  = r_6*exp(1i*k_s2*a*cos(theta_s2));
A(5,10) = - r_6*exp(-1i*k_s2*a*cos(theta_s2));

A(6,3)  = r_3*exp(1i*k_p1*a*cos(theta_p1));
A(6,4)  = - r_3*exp(-1i*k_p1*a*cos(theta_p1));
A(6,5)  = exp(1i*k_s1*a*cos(theta_s1));
A(6,6)  = exp(-1i*k_s1*a*cos(theta_s1));
A(6,7)  = - r_5*exp(1i*k_p2*a*cos(theta_p2));
A(6,8)  = r_5*exp(- 1i*k_p2*a*cos(theta_p2));
A(6,9)  = - exp(1i*k_s2*a*cos(theta_s2));
A(6,10) = - exp(-1i*k_s2*a*cos(theta_s2));

A(7,3)  = - (lambda_1 + r_3^2*(lambda_1 + 2*mu_1))*exp(1i*k_p1*a*cos(theta_p1));
A(7,4)  = - (lambda_1 + r_3^2*(lambda_1 + 2*mu_1))*exp(-1i*k_p1*a*cos(theta_p1));
A(7,5)  = - 2*mu_1*r_4*exp(1i*k_s1*a*cos(theta_s1));
A(7,6)  = 2*mu_1*r_4*exp(-1i*k_s1*a*cos(theta_s1));
A(7,7)  = (lambda_2 + r_5^2*(lambda_2 + 2*mu_2))*exp(1i*k_p2*a*cos(theta_p2));
A(7,8)  = (lambda_2 + r_5^2*(lambda_2 + 2*mu_2))*exp(-1i*k_p2*a*cos(theta_p2));
A(7,9)  = 2*mu_2*r_6*exp(1i*k_s2*a*cos(theta_s2));
A(7,10) = - 2*mu_2*r_6*exp(-1i*k_s2*a*cos(theta_s2));

A(8,3)  = - 2*mu_1*r_3*exp(1i*k_p1*a*cos(theta_p1));
A(8,4)  = 2*mu_1*r_3*exp(-1i*k_p1*a*cos(theta_p1));
A(8,5)  = mu_1*(r_4^2 - 1)*exp(1i*k_s1*a*cos(theta_s1));
A(8,6)  = mu_1*(r_4^2 - 1)*exp(-1i*k_s1*a*cos(theta_s1));
A(8,7)  = 2*mu_2*r_5*exp(1i*k_p2*a*cos(theta_p2));
A(8,8)  = - 2*mu_2*r_5*exp(-1i*k_p2*a*cos(theta_p2));
A(8,9)  = - mu_2*(r_6^2 - 1)*exp(1i*k_s2*a*cos(theta_s2));
A(8,10) = - mu_2*(r_6^2 - 1)*exp(-1i*k_s2*a*cos(theta_s2));

A(9,7)  = exp(1i*k_p2*h*cos(theta_p2));
A(9,8)  = exp(-1i*k_p2*h*cos(theta_p2));
A(9,9)  = - r_6*exp(1i*k_s2*h*cos(theta_s2));
A(9,10) = r_6*exp(-1i*k_s2*h*cos(theta_s2));
A(9,11) = - exp(1i*k_p0*h*cos(theta_p));
A(9,12) = r_2*exp(1i*k_s0*h*cos(theta_s));

A(10,7)  = r_5*exp(1i*k_p2*h*cos(theta_p2));
A(10,8)  = - r_5*exp(-1i*k_p2*h*cos(theta_p2));
A(10,9)  = exp(1i*k_s2*h*cos(theta_s2));
A(10,10) = exp(-1i*k_s2*h*cos(theta_s2));
A(10,11) = - r_1*exp(1i*k_p0*h*cos(theta_p));
A(10,12) = - exp(1i*k_s0*h*cos(theta_s));

A(11,7)  = - (lambda_2 + r_5^2*(lambda_2 + 2*mu_2))*exp(1i*k_p2*h*cos(theta_p2));
A(11,8)  = - (lambda_2 + r_5^2*(lambda_2 + 2*mu_2))*exp(-1i*k_p2*h*cos(theta_p2));
A(11,9)  = - 2*mu_2*r_6*cos(theta_s2)*exp(1i*k_s2*h*cos(theta_s2));
A(11,10) = 2*mu_2*r_6*exp(-1i*k_s2*h*cos(theta_s2));
A(11,11) = (lambda_0 + r_1^2*(lambda_0 + 2*mu_0))*exp(1i*k_p0*h*cos(theta_p));
A(11,12) = 2*mu_0*r_2*exp(1i*k_s0*h*cos(theta_s));

A(12,7)  = - 2*mu_2*r_5*exp(1i*k_p2*h*cos(theta_p2));
A(12,8)  = 2*mu_2*r_5*exp(-1i*k_p2*h*cos(theta_p2));
A(12,9)  = mu_2*(r_6^2 - 1)*exp(1i*k_s2*h*cos(theta_s2));
A(12,10) = mu_2*(r_6^2 - 1)*exp(-1i*k_s2*h*cos(theta_s2));
A(12,11) = 2*mu_0*r_1*exp(1i*k_p0*h*cos(theta_p));
A(12,12) = - mu_0*(r_2^2 - 1)*exp(1i*k_s0*h*cos(theta_s));

b(1,1) = - 1;
b(2,1) = - r_1;
b(3,1) = (lambda_0 + r_1^2*(lambda_0 + 2*mu_0));
b(4,1) = 2*mu_0*r_1;

X(:)    = linsolve(A,b);
R_p   = X(1);
R_s   = X(2).*c_p0./c_s0;
T_p1  = X(3).*c_p0./c_p1;
R_p1  = X(4).*c_p0./c_p1;
T_s1  = X(5).*c_p0./c_s1;
R_s1  = X(6).*c_p0./c_s1;
T_p2  = X(7).*c_p0./c_p2;
R_p2  = X(8).*c_p0./c_p2;
T_s2  = X(9).*c_p0./c_s2;
R_s2  = X(10).*c_p0./c_s2;
T_p   = X(11);
T_s   = X(12).*c_p0./c_s0;

J_pr = abs(R_p).^2;
J_sr = abs(R_s).^2*c_s0.*cos(theta_s)./(c_p0*cos(theta_p));
J_pt = abs(T_p).^2;
J_st = abs(T_s).^2*c_s0.*cos(theta_s)./(c_p0*cos(theta_p));