function [State_particles_new] = evolution_reaction_process(gas,State_particles,dt,...
    radiation,spark_ignition,reduced_chemistry)

Np = size(State_particles,2);

% this part is for the detailed chemistry calculation
for iNp = 1 : Np
    Psi_iNp = State_particles(:,iNp);
    hmass_iNp = Psi_iNp(1,1);
    P_iNp = Psi_iNp(2,1);
    Y_iNp = Psi_iNp(3:end,1);
    set(gas, 'Enthalpy', hmass_iNp, 'Pressure', P_iNp,'MassFractions',Y_iNp);
    Psi_0=[temperature(gas);
        P_iNp;
        Y_iNp];
    if true(reduced_chemistry.condition)
        options = odeset('Mass',reduced_chemistry.Ms,'RelTol',1.e-7,'AbsTol',1.e-9);
    else
        options = odeset('RelTol',1.e-8,'AbsTol',1.e-10);
    end
    out = ode15s(@conhp,[0 dt],Psi_0,options,gas,radiation,spark_ignition);
    %
    set(gas, 'Temperature', out.y(1,end), 'Pressure', out.y(2,end),...
        'MassFractions',out.y(3:end,end));
    State_particles_new(:,iNp) = [enthalpy_mass(gas);
        pressure(gas);
        massFractions(gas) ];
end


end

function F = conhp(t,Psi,gas,radiation,spark_ignition)

% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, Psi(3:end,1), 'nonorm');
set(gas, 'Temperature', Psi(1,1), 'Pressure', Psi(2,1));

% get molecular weights, [mw]=[kg/kmol];
mw = molecularWeights(gas);

% energy equation
omega_dot = netProdRates(gas);

%%%%%%%%%%%%%%%%%
T_dot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
    * omega_dot / (density(gas)*cp_mass(gas));

if true(radiation.condition)
    X = moleFractions(gas);
    for i = 1 : size(radiation.radiating_species,2)
        X_radiating_species(i) = X(speciesIndex(gas,radiation.radiating_species{i}));
    end
    [q_rad_volumetric] = get_radiative_heat_loss(temperature(gas),X_radiating_species,...
        pressure(gas),radiation.radiating_species);
    T_dot = T_dot - q_rad_volumetric/ (density(gas)*cp_mass(gas));
end

if true(spark_ignition.condition)
    if spark_ignition.time < spark_ignition.spark_duration_time
        T_dot = T_dot + spark_ignition.qs_volumetric / (density(gas)*cp_mass(gas));
    end
end


% pressure equation
P_dot = 0;

% species equations
% get the inverse of density, [rho]=kg/m3
rrho = 1.0/density(gas);
dYdt = rrho*(mw.*omega_dot);

% set up column vector for dydt
F = [ T_dot;
    P_dot;
    dYdt ];

end
