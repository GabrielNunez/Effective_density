function theta_s = Reflection_angle_s_waves(lambda, mu, rho_x, rho_z, theta_i, delta, tol)

% if mu < 1e0
%     mu = vpa(mu);
% end

c_px = sqrt((lambda + 2*mu)/rho_x);

n = 0;
ig    = 1/c_px;
err   = 1;

[c_1, ~,  ~] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda, mu, rho_x, rho_z, ig);
c_s  = c_1(2);

[c_1, ~,  ~] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda, mu, rho_x, rho_z, theta_i);
c_p  = c_1(1);

while err > tol && n < 10

    f_0 = sin(ig)/c_s - sin(theta_i)/c_p;

    [c_1, ~,  ~] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda, mu, rho_x, rho_z, ig - delta);
    c_s_d  = c_1(2);
    
    f_0_d = sin(ig - delta)/c_s_d - sin(theta_i)/c_p;

    f_p0 = (f_0_d - f_0)/delta;

    ig = ig + f_0/f_p0;

    [c_1, ~,  ~] = Anisotropic_density_phase_velocities_polarisation_vectors(lambda, mu, rho_x, rho_z, ig);
    c_s  = c_1(2);

    err = abs((sin(ig)/c_s - sin(theta_i)/c_p)/(sin(theta_i)/c_p));

    n = n + 1; 
end

theta_s = ig;