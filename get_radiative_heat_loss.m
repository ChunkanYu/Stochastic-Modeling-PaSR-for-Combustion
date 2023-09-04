function [q_rad_volumetric] = get_radiative_heat_loss(T,X,p_total,radiation_species)
% [q_rad_volumetric] = W/m^3
sigma = 5.669e-8;

k = 0;

for i = 1 : size(radiation_species,2)
    a_pi = Planck_mean_absorption_coefficient(T,radiation_species{i});
    p_i = X(i) * p_total *1e-5; % transfer to [bar]
    k = k + a_pi * p_i;
end

q_rad_volumetric = 4 * k * sigma * (T^4 - 300^4);

end





function [a_p] = Planck_mean_absorption_coefficient(T,species)
% unit 1/(m*bar)

switch species
    case 'NO'
        if (300<=T) &&(T<=785)
            a_p = 3.05727e1 + (-3.71207e-1)*T + (1.76286e-3)*T^2 - ...
                (4.14129e-6)*T^3 + (5.24120e-9)*T^4 - (3.45099e-12)*T^5 + ...
                (9.32118e-16)*T^6;
        else
            a_p = 3.28579e0 + (4.89953e-3)*T - (1.91582e-5)*T^2 + ...
                (1.98303e-8)*T^3 - (9.75832e-12)*T^4 + (2.37376e-15)*T^5 - ...
                (2.30104e-19)*T^6;
        end
    % 
    case 'N2O'
        if (300<=T) &&(T<=820)
            a_p = 9.34925e1 + (-1.04332e0)*T + (4.84016e-3)*T^2 - ...
                (1.03506e-5)*T^3 + (1.15054e-8)*T^4 - (6.58670e-12)*T^5 + ...
                (1.55603e-15)*T^6;
        else
            a_p = 2.69054e1 + (1.24925e-1)*T - (3.70499e-4)*T^2 + ...
                (3.73795e-7)*T^3 - (1.83887e-10)*T^4 + (4.49341e-14)*T^5 - ...
                (4.37851e-18)*T^6;
        end
    %
    case 'H2O'
        a_p = (-1.38759e0) + (5.07308e3)/T - (2.45950e6)/T^2 + ...
            (5.80188e9)/T^3 - (2.30409e12)/T^4 + (3.24628e14)/T^5;
    %
    case 'NH3'
        if (300<=T) &&(T<=1100)
            a_p = (1.14178e1) - (3.51460e4)/T + (3.29078e7)/T^2 - ...
            (6.99472e9)/T^3 + (2.23403e11)/T^4 + (6.18633e13)/T^5;
        else
            a_p = (-5.04091e-2) - (8.49774e2)/T + (9.11833e6)/T^2 - ...
            (2.90195e10)/T^3 + (3.60444e13)/T^4 - (1.27775e16)/T^5;
        end
    %
end

end