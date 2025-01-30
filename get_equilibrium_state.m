function [Psi_eq, T_eq]=get_equilibrium_state(gas)
% calculate the equilibrium state while keeping the enthalpy and pressure
% constant
%
% output: Psi_eq = (h,p,wi);
%                   h: specific enthalpy [J/kg]
%                   p: pressure [Pa]
%                   wi: mass fraction [-] 


equilibrate(gas,'HP'); 

hmass_eq = enthalpy_mass(gas);

T_eq = temperature(gas);

P_eq= pressure(gas);

Y_eq = massFractions(gas);

X_eq = moleFractions(gas);

Psi_eq = [hmass_eq;
          P_eq;
          Y_eq];

      
      
end

