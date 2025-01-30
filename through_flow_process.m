function [State_particles,RefVar_particles]=through_flow_process(flame_type,...
    State_particles, State_unburnt, N_replace,mixing_model_type,...
    RefVar_particles)


Np = size(State_particles,2);


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
        State_particles(:,ip_particle_replace(i)) = State_unburnt;
        end
    case 'nonpremixed'
        if N_replace == 1
          if (rand(1)<0.5)
              State_particles(:,ip_particle_replace(1)) = State_unburnt.oxidizer;
          else
              State_particles(:,ip_particle_replace(1)) = State_unburnt.fuel;
          end
        else
        N_replace_oxidizer = round(0.5*N_replace);
        for i=1:N_replace_oxidizer
            State_particles(:,ip_particle_replace(i)) = State_unburnt.oxidizer;
        end
        for i = N_replace_oxidizer + 1: N_replace
            State_particles(:,ip_particle_replace(i)) = State_unburnt.fuel;            
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