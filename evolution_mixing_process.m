function [State_particles,RefVar_particles]=evolution_mixing_process(mixing_model_type,...
    State_particles, RefVar_particles,mixing_model_parameter,omega_turb,dt)

switch mixing_model_type
    case 'IEM'
        [State_particles]=evolution_mixing_process_IEM(State_particles,...
            mixing_model_parameter,omega_turb,dt);
    case 'MCM'
        [State_particles]=evolution_mixing_process_MCM(State_particles,...
            mixing_model_parameter,omega_turb,dt);
    case 'MMC_MCM_Sundaram'
        [State_particles,RefVar_particles]=evolution_mixing_process_MMC_MCM_Sundaram(State_particles,...
            RefVar_particles,mixing_model_parameter,omega_turb,dt);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IEM
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MCM
function [Psi_particles]=evolution_mixing_process_MCM(Psi_particles,...
    mixing_model_parameter,omega_turb,dt)

C_phi = mixing_model_parameter.C_phi;

Np = size(Psi_particles,2);

mix_num = 3 * Np * C_phi * omega_turb * dt;
for i = 1 : mix_num
    % select 2 particles randomly from Np - particles
    ip_particle_select = randperm( Np,2 );
    % determine whether the particles should be mixed with each other
    if (rand(1) < 0.5)
        % mixing extent
        MCM_extent = rand(1);
        % mean value of two particles
        mean_Psi_particles_selected = 0.5 * sum(Psi_particles(:,ip_particle_select),2);
        % evolution of particle thermo-kinetic state
        Psi_particles(:,ip_particle_select(1)) = (1 - MCM_extent) * Psi_particles(:,ip_particle_select(1)) + ...
            MCM_extent * mean_Psi_particles_selected;

        Psi_particles(:,ip_particle_select(2)) = (1 - MCM_extent) * Psi_particles(:,ip_particle_select(2)) + ...
            MCM_extent * mean_Psi_particles_selected;
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MMC - MCM-Sundaram
function [Psi_particles,RefVar_particles]=evolution_mixing_process_MMC_MCM_Sundaram(Psi_particles,...
    RefVar_particles,mixing_model_parameter,omega_turb,dt)

% evolution of reference variable
% effective drift velocity
u = mixing_model_parameter.L / mixing_model_parameter.residence_time;

% effective diffusion 
b = mixing_model_parameter.lambda * u;

Np = size(RefVar_particles,2);

% Ornstein-Uhlenbeck process
for i = 1 : Np
    RefVar_particles(i) = RefVar_particles(i) + u * dt + b* sqrt(dt) * randn(1);
end

condition = 'false';

while (condition=='false')
% make sure that the RefVar greater than 0 by reflecting at x=0
RefVar_particles = abs(RefVar_particles);

% make sure that the RefVar less than L by reflecting at x=L
[ind_outside_domain] = find(RefVar_particles > mixing_model_parameter.L);
RefVar_particles(ind_outside_domain) = 2*mixing_model_parameter.L - ...
    RefVar_particles(ind_outside_domain);

[ind] = find(RefVar_particles<0 | RefVar_particles > mixing_model_parameter.L);
condition = isempty(ind);

end


% mixing based on reference variable
[~,ind] = sort(RefVar_particles);

C_phi = mixing_model_parameter.C_phi;

Np = size(Psi_particles,2);

mix_num = 3 * Np * C_phi * omega_turb * dt;

for i = 1 : mix_num
    % select 1 particle randomly from Np - particles
    ip_particle_select(1) = randperm( Np,1 );
    % find the next particle which is close 
    [ind_ip]=find(ind == ip_particle_select(1));
    if ind_ip == 1
        ip_particle_select(2) = ind(ind_ip+1);
    else
    ip_particle_select(2) = ind(ind_ip-1);
    end
    if (rand(1) < 0.5)
        % mixing extent
        MCM_extent = rand(1);
        % mean value of two particles
        mean_Psi_particles_selected = 0.5 * sum(Psi_particles(:,ip_particle_select),2);
        % evolution of particle thermo-kinetic state
        Psi_particles(:,ip_particle_select(1)) = (1 - MCM_extent) * Psi_particles(:,ip_particle_select(1)) + ...
            MCM_extent * mean_Psi_particles_selected;

        Psi_particles(:,ip_particle_select(2)) = (1 - MCM_extent) * Psi_particles(:,ip_particle_select(2)) + ...
            MCM_extent * mean_Psi_particles_selected;
    end
end


end
