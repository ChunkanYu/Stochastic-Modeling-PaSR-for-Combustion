function [w_species_SS,Temperature_SS] = homogeneous_reactor(Psi_unburnt,gas,reactor_time)


reactor_system='isochor';
chemistry='detailed_chemistry'; 

source_term.condition = false;
source_term.E_spark = 1.1e9;
source_term.spark_duration = 1e-3;


mw = molecularWeights(gas);
set(gas,'Enthalpy',Psi_unburnt(1,1),'Pressure',Psi_unburnt(2,1),...
    'MassFractions',Psi_unburnt(3:end,1));

nsp = nSpecies(gas);


y0 = [temperature(gas)
    massFractions(gas)];


tel = [0 reactor_time];

Ms=eye(nsp+1,nsp+1);

switch chemistry
    case 'detailed chemistry'
    case 'GQL_chemistry'
        Ms = importdata('GQL_Ms.mat');
%         Ms = Ms_GQL;
    case 'QSSA_chemistry'
        for i = 1 : size(QSS_species,2)
            iQSS_species = speciesIndex(gas,QSS_species{i});
            if iQSS_species ~= 0
            Ms(iQSS_species+1,iQSS_species+1) = 0;
            else
                fprintf([QSS_species{i},'-QSS is not found in this mechanism\n']);
            end
        end
end



options = odeset('Mass',Ms(:,:,1),'RelTol',1.e-8,'AbsTol',1.e-10);

out = ode15s(@ode_rhs,tel,y0,options,gas,mw,reactor_system,source_term);

Temperature_SS = out.y(1,end);
w_species_SS = out.y(2:end,end);
pressure_SS = 


end



function dydt = ode_rhs(t, y, gas, mw,reactor_system,source_term)


% Set the state of the gas, based on the current solution vector.
setMassFractions(gas,y(2:end),'nonorm')

nsp = nSpecies(gas);


if strcmp(reactor_system,'isochor')
    set(gas, 'T', y(1), 'Rho', density(gas));
    % energy equation
    wdot = netProdRates(gas);
    tdot = - temperature(gas) * gasconstant * (enthalpies_RT(gas) - ones(nsp,1))' ...
        * wdot / (density(gas)*cv_mass(gas));
    % species equation
    rrho = 1.0/density(gas);
    for i = 1:nsp
        dwdt(i,1) = rrho*mw(i)*wdot(i);
    end
end

if strcmp(reactor_system,'isobar') 
    set(gas, 'T', y(1), 'P', pressure(gas));
    % energy equation
    wdot = netProdRates(gas);
    
    tdot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
        * wdot / (density(gas)*cp_mass(gas));

    if true(source_term.condition)
        if t<=source_term.spark_duration
            tdot = tdot + source_term.E_spark/(density(gas)*cp_mass(gas));
        end
    end


    % species equation
    rrho = 1.0/density(gas);
    for i = 1:nsp
        dwdt(i,1) = rrho*mw(i)*wdot(i);
    end

    

end

% set up column vector for dydt
dydt = [ tdot;
    dwdt];
end
