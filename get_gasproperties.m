function [gas] = get_gasproperties()
%GET_GASPROPERTIES Look up gas properties using python library coolprop.
%Necessit coolprop installed
%   Detailed explanation goes here
%% perfect gas caracteristics
gas.composition = 'methane[0.9318]&ethane[0.0435]&nitrogen[0.0096]&CO2[0.0078]&propane[0.0073]';
gas.T_st = 273.15; %standard conditions temperature (K)
gas.T_avg = 288.15; %avg. gound temperature = 15°C (K)
gas.p_st = 1.01325e5; %standard conditions pressure (Pa)
gas.rho_st = py.CoolProp.CoolProp.PropsSI('D','P',gas.p_st,'T',gas.T_st,gas.composition); %density in kg/m^3
gas.rho_air = py.CoolProp.CoolProp.PropsSI('D','P',gas.p_st,'T',gas.T_st,'air');
gas.d_rel = gas.rho_st/gas.rho_air;
gas.HHV = 41.328e6; %Higher Heating Value in J/m^3 st
gas.mu = py.CoolProp.CoolProp.PropsSI('viscosity','P',gas.p_st,'T',gas.T_st,gas.composition); %dynamic viscosity in Pa*s
gas.Molar_mass = py.CoolProp.CoolProp.PropsSI('M','P',gas.p_st,'T',gas.T_st,gas.composition); %density in kg/m^3
gas.R_specific = 8.31432/gas.Molar_mass/1000; %Specific gas constant in kJ/kg/K
gas.z_st = py.CoolProp.CoolProp.PropsSI('Z','P',gas.p_st,'T',gas.T_st,gas.composition); %compressibility factor
gas.H = py.CoolProp.CoolProp.PropsSI('H','P',gas.p_st,'T',gas.T_st,gas.composition)/1000; %specific enthalpy in kJ/kg
gas.Cp = py.CoolProp.CoolProp.PropsSI('CP0MASS','P',gas.p_st,'T',gas.T_st,gas.composition)/1e3; %Ideal gas mass specific constant pressure specific heat in kJ/kg/K
gas.v = py.CoolProp.CoolProp.PropsSI('viscosity','P',gas.p_st,'T',gas.T_st,gas.composition)*1e6; %Viscosity in uPa*s
gas.cond = py.CoolProp.CoolProp.PropsSI('conductivity','P',gas.p_st,'T',gas.T_st,gas.composition)*1e3; %Thermal conductivity in mW/m/K

gas.z_avg = 1;

gas.MW_to_m3 = 1e6./(gas.rho_st*gas.HHV);
gas.m3_to_MW = 1/gas.MW_to_m3;


%c = sqrt(p_ref/rho_ref); %speed of sound, assumed constant
end

