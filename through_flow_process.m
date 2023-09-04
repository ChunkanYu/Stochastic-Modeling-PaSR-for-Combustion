function [Psi_particles,RefVar_particles]=through_flow_process(flame_type,...
    Psi_particles, Psi_unburnt, N_replace,mixing_model_type,...
    RefVar_particles)


Np = size(Psi_particles,2);


% pointer to the particles needed to be replaced
if strcmp(mixing_model_type,'MMC_MCM_Sundaram')
    [~,ind_sort]=sort(RefVar_particles(1,:),'descend');
    ip_particle_replace = ind_sort(1,1:N_replace);
else
    ip_particle_replace = randperm(Np,N_replace);
end


switch flame_type
    case 'premixed'
        for i=1:N_replace
        Psi_particles(:,ip_particle_replace(i)) = Psi_unburnt;
        end
    case 'nonpremixed'
        if N_replace == 1
          if (rand(1)<0.5)
              Psi_particles(:,ip_particle_replace(1)) = Psi_unburnt.oxidizer;
          else
              Psi_particles(:,ip_particle_replace(1)) = Psi_unburnt.fuel;
          end
        else
        N_replace_oxidizer = round(0.5*N_replace);
        for i=1:N_replace_oxidizer
            Psi_particles(:,ip_particle_replace(i)) = Psi_unburnt.oxidizer;
        end
        for i = N_replace_oxidizer + 1: N_replace
            Psi_particles(:,ip_particle_replace(i)) = Psi_unburnt.fuel;            
        end
        end
end

for i=1:N_replace
    switch mixing_model_type
        case 'MMC_IEM_Varna'
            RefVar_particles(1,ip_particle_replace(i)) = 0; % 0 for the unburnt state
        case 'MMC_MCM_Sundaram'
            RefVar_particles(1,ip_particle_replace(i)) = 0;
    end
end









end