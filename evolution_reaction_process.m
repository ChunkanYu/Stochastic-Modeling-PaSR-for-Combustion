function [Psi_particles_new,Temperature_particles] = evolution_reaction_process(gas,Psi_particles,dt,...
    radiation)


Np = size(Psi_particles,2);

for i=1:Np
    Psi_i = Psi_particles(:,i);
    hmass_i = Psi_i(1,1);
    P_i = Psi_i(2,1);
    Y_i = Psi_i(3:end,1);
    set(gas, 'Enthalpy', hmass_i, 'Pressure', P_i,'MassFractions',Y_i);
    Psi_0=[temperature(gas);
        P_i;
        Y_i];
    options = odeset('RelTol',1.e-8,'AbsTol',1.e-10);
    out = ode15s(@conhp,[0 dt],Psi_0,options,gas,radiation);
    
    %
    set(gas, 'Temperature', out.y(1,end), 'Pressure', out.y(2,end),...
        'MassFractions',out.y(3:end,end));
    Psi_particles_new(:,i) = [enthalpy_mass(gas);
        pressure(gas);
        massFractions(gas) ];
    Temperature_particles(i) = temperature(gas);
end



end

function F = conhp(t,Psi,gas,radiation)

% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, Psi(3:end,1), 'nonorm');
set(gas, 'Temperature', Psi(1,1), 'Pressure', Psi(2,1));

% energy equation
omega_dot = netProdRates(gas);
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


% pressure equation
P_dot = 0;

% species equations
dYdt = reaction_rate(gas);

% set up column vector for dydt
F = [ T_dot;
    P_dot;
    dYdt ];

end

function dYdt = reaction_rate(gas)
% get the species equation which fulfills:
% dY/dt = F;

% get molecular weights, [mw]=[kg/kmol];
mw = molecularWeights(gas);

% get the species production rates, [wdot]=[kmol/m3 s]
omega_dot = netProdRates(gas);

% get the inverse of density, [rho]=kg/m3
rrho = 1.0/density(gas);

% get the species equation, [1/s]
dYdt = rrho*(mw.*omega_dot);
end