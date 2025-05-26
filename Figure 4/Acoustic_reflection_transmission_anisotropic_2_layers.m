function [R, T, R_1, T_1, R_2, T_2, I_r, I_t, I_r1, I_t1, I_r2, I_t2] = Acoustic_reflection_transmission_anisotropic_2_layers(freq, rho_0x, rho_0z, rho_1, rho_2, K_0, K_1, K_2, h, varphi, theta_i)

omega = 2*pi*freq;

c_0 = sqrt(K_0*(sin(theta_i).^2/rho_0x + cos(theta_i).^2/rho_0z));
c_1 = sqrt(K_1/rho_1);
c_2 = sqrt(K_2/rho_2);

theta_t1 = asin(c_1./c_0.*sin(theta_i));
theta_t2 = asin(c_2./c_0.*sin(theta_i));

z_0 = rho_0z*c_0./cos(theta_i);
z_1 = rho_1*c_1./cos(theta_t1);
z_2 = rho_2*c_2./cos(theta_t2);

k_0   = omega./c_0;
k_1   = omega./c_1;
k_2   = omega./c_2;

k_0x = k_0.*cos(theta_i);
k_1x = k_1.*cos(theta_t1);
k_2x = k_2.*cos(theta_t2);

a = varphi*h;
b = (1 - varphi)*h;

D = 2*z_0.*z_1.*z_2.*cos(k_1x*a).*cos(k_2x*b) - z_0.*(z_2.^2 + z_1.^2).*sin(k_1x*a).*sin(k_2x*b) - 1i*z_1.*(z_2.^2 + z_0.^2).*cos(k_1x*a).*sin(k_2x*b) - 1i*z_2.*(z_1.^2 + z_0.^2).*sin(k_1x*a).*cos(k_2x*b);

R   = (z_0.*(z_2.^2 - z_1.^2).*sin(k_1x*a).*sin(k_2x*b) - 1i*z_1.*(z_2.^2 - z_0.^2).*cos(k_1x*a).*sin(k_2x*b) - 1i*z_2.*(z_1.^2 - z_0.^2).*sin(k_1x*a).*cos(k_2x*b))./D;
T_1 = z_1.*exp(+1i*k_1*a).*(z_2.*(z_0 + z_1).*cos(k_2.*b) + 1i*(z_2.^2 + z_0.*z_1).*sin(k_2*b))./D;
R_1 = z_1.*exp(-1i*k_1*a).*(z_2.*(z_0 - z_1).*cos(k_2.*b) + 1i*(z_2.^2 - z_0.*z_1).*sin(k_2*b))./D;
T_2 = z_1.*z_2*exp(+1i*k_2*h).*(z_0 + z_2)./D;
R_2 = z_1.*z_2*exp(-1i*k_2*h).*(z_0 - z_2)./D;
T   = 2.*z_0.*z_1.*z_2.*exp(-1i*k_0x*h)./D;

I_r  = real(abs(R).^2);
I_t  = real(abs(T).^2);
I_r1 = real(abs(R_1).^2.*z_0./z_1);
I_t1 = real(abs(T_1).^2.*z_0./z_1);
I_r2 = real(abs(R_2).^2.*z_0./z_2);
I_t2 = real(abs(T_2).^2.*z_0./z_2);