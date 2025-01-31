clear all;clc; close all;
% function [ignition] = main_Stochastic_PaSR
Np=50; % number of particle % all 600 Np
% nReac = 0;
T_unburnt = 300.0; % temperature of unburnt gas
p0 = 1.e5; % pressure of the system
Phi = 1.0; % fuel/air equivalent ratio
alpha_h = 0.3; % hydrogen content in NH3-H2-mixture

C_phi = 2.0; % mixing model parameter
% the following for MMC-MCM-Sundaram
L = 0.1; lambda = 0.1;

dt=2e-5; % all 10/6
residence_time = 1.0e-3;

flame_type = 'premixed';
% mixing_model_type = 'MMC_MCM_Sundaram';
mixing_model_type = 'MCM';
omega_turb = 3.e3; % 3e3 with 1ms no extinction

gas = Solution('../mechanism_H2_Air/Warnatz.cti');
nsp = nSpecies(gas);

% set-up block for radiation effect
radiation.condition = false;
radiation.radiating_species={'NO','N2O','H2O','NH3'};

% set-up block for spark ignition process
spark_ignition.condition = false;
spark_ignition.qs_volumetric = 6e8;
spark_ignition.spark_duration_time = 1e-3;


% set-up block for DAE reduced chemistry (e.g. QSSA, GQL)
reduced_chemistry.condition = false;
reduced_chemistry.chemistry='QSSA_chemistry'; %
QSS_species={'OH'};

% set-up block for CRN / pre-chamber
prechamber.condition = false;
prechamber.reactor_time = 1e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ih2o = speciesIndex(gas,'H2O');
io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
in2o = speciesIndex(gas,'N2O');
ico = speciesIndex(gas,'CO');

mw = molecularWeights(gas);

[Psi_unburnt,gas]=get_unburnt_state(gas,T_unburnt,p0,alpha_h,Phi,flame_type);

if true(prechamber.condition)
  [w_species_SS,Temperature_SS] = homogeneous_reactor(Psi_unburnt,gas,prechamber.reactor_time);
end

[Psi_eq, T_eq]=get_equilibrium_state(gas);

Psi_particles = repmat(Psi_eq,1,Np);
% Psi_particles = repmat(0.5*Psi_eq + 0.5*Psi_unburnt,1,Np);

RefVar_particles = zeros(1,Np);


switch mixing_model_type
    case 'MMC_MCM_Sundaram'
        RefVar_particles = 0 + (L - 0)*rand(1,Np);
        mixing_model_parameter.lambda = lambda;
        mixing_model_parameter.L = L;
        mixing_model_parameter.residence_time = residence_time;
end

mixing_model_parameter.C_phi = C_phi;

% here determine the number of particles to be replaced
N_replace = floor( Np * dt / residence_time );

if true(spark_ignition.condition)
    State_particles = repmat(Psi_unburnt,1,Np);
else  
    State_particles = Psi_particles;
end

State_unburnt = Psi_unburnt;

% here prepare for the reduced chemistry
%       QSSA
if true(reduced_chemistry.condition)
    switch reduced_chemistry.chemistry
        case 'GQL_chemistry'
            reduced_chemistry.Ms = importdata('GQL_Ms.mat');
        case 'QSSA_chemistry'
            reduced_chemistry.Ms =eye(nsp+2,nsp+2);
            for i = 1 : size(QSS_species,2)
                iQSS_species = speciesIndex(gas,QSS_species{i});
                if iQSS_species ~= 0
                    reduced_chemistry.Ms(iQSS_species+2,iQSS_species+2) = 0;
                else
                    fprintf([QSS_species{i},'-QSS is not found in this mechanism\n']);
                end
            end
    end
end




for n = 1: 1e5
    % through flow
    [State_particles,RefVar_particles]=through_flow_process(flame_type,...
        State_particles, State_unburnt, N_replace,mixing_model_type,...
        RefVar_particles);
    % if mod(n,2) == 0
    %     [Psi_particles,RefVar_particles]=through_flow_process(flame_type,...
    %     Psi_particles, Psi_unburnt, 1, mixing_model_type,...
    %     RefVar_particles);
    % end
    % mixing process
    % omega_turb = (2e4-1e3).*rand(1,1) + 1e3;
    [State_particles,RefVar_particles]=evolution_mixing_process(mixing_model_type,...
        State_particles, RefVar_particles,mixing_model_parameter,omega_turb,dt);
    % reaction process
    if true(spark_ignition.condition)
        spark_ignition.time = (n-1) * dt;
    end
    [State_particles] = evolution_reaction_process(gas,State_particles,dt,...
        radiation,spark_ignition,reduced_chemistry);

    %

    particles_hpwi(:,:,n) = State_particles;
    for iNp = 1 : Np
        set(gas, 'Enthalpy', particles_hpwi(1,iNp,n), 'Pressure', particles_hpwi(2,iNp,n),...
            'MassFractions',particles_hpwi(3:end,iNp,n));
        particles_temperature(iNp) = temperature(gas);
    end
    % plot

    for k=1:Np
        YYY = particles_hpwi(3:end,k,n)./mw;
        mole_fraction(:,k) = (1/sum(YYY))*YYY;
    end
    mean_TT(n) = sum(particles_temperature)/Np;
    mean_all(:,n)=(1/Np)*sum(mole_fraction,2);

    plot(mean_TT,'ro'); hold on; pause(0.01);
    
    % if spark_ignition.time > spark_ignition.spark_duration_time
    %     if mean_TT(n) > 1200
    %     ignition = 1;
    %     else
    %         ignition = 0;
    %     end
    %     break
    % end
end


% mean_NO = sum(mean_all(ino,300:450))/151;

% close all;
