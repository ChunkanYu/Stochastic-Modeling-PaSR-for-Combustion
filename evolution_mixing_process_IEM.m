function [Psi_particles]=evolution_mixing_process_IEM(Psi_particles,...
    mixing_model_parameter,omega_turb,dt)
% Psi_particle: nsp x Np matrix
%    |     |         |
% [Psi_1 Psi_2 *** Psi_Np]
%    |     |         |
% dt: time step

% Block for the Interaction by Exchange with the Mean (IEM) model


C_phi = mixing_model_parameter.C_phi;


Np = size(Psi_particles,2);
Psi_mean = sum(Psi_particles,2)/Np;
for i = 1 : Np
    Psi_particles(:,i) = Psi_mean + ...
        (Psi_particles(:,i)-Psi_mean)*exp(-0.5*C_phi*omega_turb*dt);
end



end