clear all;clc; close all;
% function [mean_NO] = main_Stochastic_PaSR(nReac)
Np=400; % number of particle % all 600 Np
% nReac = 0;
T_unburnt = 300.0; % temperature of unburnt gas
p0 = 40.e5; % pressure of the system
Phi = 1.0; % fuel/air equivalent ratio
alpha_h = 0.3; % hydrogen content in NH3-H2-mixture

C_phi = 2.0; % mixing model parameter
% the following for MMC-IEM-Varna
C_RefVar = 2.0; b0 = 0.71; C_min = 8.0;
% the following for MMC-MCM-Sundaram
L = 0.1; lambda = 0.1;

dt=1.5*(10/6)*1e-6; % all 10/6
residence_time = 1.0e-3; 

flame_type = 'premixed';
% mixing_model_type = 'MMC_MCM_Sundaram';
mixing_model_type = 'IEM';
omega_turb = 1.e5; % 2e3 with 1ms extinction

radiation.condition = true;
radiation.radiating_species={'NO','N2O','H2O','NH3'};


gas = Solution('../Mech_NH3_H2/shrestha2021_noncarbon_chem.cti');
% gas = Solution('../Mech_H2/Li_fuel_2019_only_h2.cti');
% gas = Solution('../Mech_NH3_H2/Nakamura_NH3_H2.cti');
% gas = Solution('../Mech_NH3_H2/Stagni_NH3_H2.cti');
% gas = Solution('../Mech_NH3_H2/Gotama_NH3_H2.cti');
% 
% if nReac ~=0
% setMultiplier(gas,nReac,2);
% end


ih2o = speciesIndex(gas,'H2O');
ino = speciesIndex(gas,'NO');
ino2 = speciesIndex(gas,'NO2');
in2o = speciesIndex(gas,'N2O');
ico = speciesIndex(gas,'CO');

[Psi_unburnt,gas]=get_unburnt_state(gas,T_unburnt,p0,alpha_h,Phi,flame_type);

[Psi_eq, T_eq]=get_equilibrium_state(gas);


Psi_particles = repmat(Psi_eq,1,Np);
RefVar_particles = zeros(1,Np);

switch mixing_model_type
    case 'MMC_IEM_Varna'
        RefVar_particles = ones(1,Np);
        mixing_model_parameter.C_RefVar = C_RefVar;
        mixing_model_parameter.b0 = b0;
        mixing_model_parameter.C_min = C_min;

    case 'MMC_MCM_Sundaram'
        RefVar_particles = 0 + (L - 0)*rand(1,Np);
        mixing_model_parameter.lambda = lambda;
        mixing_model_parameter.L = L;
        mixing_model_parameter.residence_time = residence_time;
end


mixing_model_parameter.C_phi = C_phi;


mw = molecularWeights(gas);

% here determine the number of particles to be replaced
N_replace = floor( Np * dt / residence_time );

for n = 1: 1e5
    % through flow
    [Psi_particles,RefVar_particles]=through_flow_process(flame_type,...
        Psi_particles, Psi_unburnt, N_replace,mixing_model_type,...
        RefVar_particles);
    % if mod(n,2) == 0
    %     [Psi_particles,RefVar_particles]=through_flow_process(flame_type,...
    %     Psi_particles, Psi_unburnt, 1, mixing_model_type,...
    %     RefVar_particles);
    % end
    % mixing process
    switch mixing_model_type
        case 'IEM'
            [Psi_particles]=evolution_mixing_process_IEM(Psi_particles,...
                mixing_model_parameter,omega_turb,dt);
        case 'MCM'
            [Psi_particles]=evolution_mixing_process_MCM(Psi_particles,...
                mixing_model_parameter,omega_turb,dt);
        case 'MMC_IEM_Varna'
            [Psi_particles,RefVar_particles]=evolution_mixing_process_MMC_IEM_Varna(Psi_particles,...
                RefVar_particles,mixing_model_parameter,omega_turb,dt);
        case 'MMC_MCM_Sundaram'
            [Psi_particles,RefVar_particles]=evolution_mixing_process_MMC_MCM_Sundaram(Psi_particles,...
                RefVar_particles,mixing_model_parameter,omega_turb,dt);
    end
    % reaction process
    [Psi_particles,TT] = evolution_reaction_process(gas,Psi_particles,dt,radiation);
    
    % 
    particles_hpwi(:,:,n) = Psi_particles;
    % plot

    for k=1:Np
        YYY = Psi_particles(3:end,k)./mw;
        mole_fraction(:,k) = (1/sum(YYY))*YYY;
    end
    mean_TT(n) = sum(TT)/Np;
    mean_all(:,n)=(1/Np)*sum(mole_fraction,2);
    mean_NOx(n) = mean_all(ino,n)+mean_all(ino2,n);
    subplot(1,2,1)
    semilogy(n,1e6*mean_NOx(n),'ko'); hold on;
    subplot(1,2,2)
    plot((n-1)*dt,mean_TT(n),'ro'); hold on;
    pause(.01);
end


% mean_NO = sum(mean_all(ino,300:450))/151;

% close all;
