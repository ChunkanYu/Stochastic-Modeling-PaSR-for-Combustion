function [Psi_unburnt,gas]=get_unburnt_state(gas,T_unburnt,p0,...
    alpha_h,Phi,flame_type)
% calculate the equilibrium state while keeping the enthalpy and pressure
% constant
%
% output: Psi_ini = (h,p,wi);
%                   h: specific enthalpy [J/kg]
%                   p: pressure [Pa]
%                   wi: mass fraction [-]


nsp = nSpecies(gas);
X = zeros(nsp,1);

io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
ih2 = speciesIndex(gas,'H2');
inh3 = speciesIndex(gas,'NH3');
ich4 = speciesIndex(gas,'CH4');



% NH3 - based mixture
% NH3 + H2
X(io2) = (0.75 - alpha_h/4)/Phi;
X(in2) = (0.75 - alpha_h/4)/Phi*(79/21);
X(ih2) = alpha_h;
X(inh3) =  1- alpha_h;

% X(io2) = 0.75;
% X(in2) = 1.75;
% X(inh3) = 1.0;

switch flame_type
    case 'premixed'
        set(gas,'Temperature',T_unburnt,'Pressure',p0,'MoleFractions',X);
        Psi_unburnt = [enthalpy_mass(gas);
            pressure(gas);
            massFractions(gas)];
    case 'nonpremixed'
        X_oxidizer = zeros(nsp,1); X_fuel = zeros(nsp,1);
        X_oxidizer(io2) = 1; X_oxidizer(in2) = 79/21;X_fuel(ich4) = 1.0;
        set(gas,'Temperature',T_unburnt,'Pressure',p0,'MoleFractions',X_oxidizer);
        Psi_unburnt.oxidizer = [enthalpy_mass(gas);
                                pressure(gas);
                                massFractions(gas)];

        set(gas,'Temperature',T_unburnt,'Pressure',p0,'MoleFractions',X_fuel);

        Psi_unburnt.fuel = [enthalpy_mass(gas);
                            pressure(gas);
                            massFractions(gas)];

        set(gas,'Temperature',T_unburnt,'Pressure',p0,'MoleFractions',X);
end

end

