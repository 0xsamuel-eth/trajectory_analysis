clc;clear;
mp = 42.5;
m_ch4 = 8.5;
m_o2 = 34;
mm_ch4 = 16;
mm_o2 = 32;

frac_ch4 = m_ch4/mp;
frac_o2 = m_o2/mp;

mwmix = 1/((frac_ch4/mm_ch4)+(frac_o2/mm_o2)); %kg/kmol
ru = 8314.5;

R = ru/mwmix;


massflow = mp/7;

Ptotal = 25* 101325;
gamma = 1.4;
ttotal = 3533; % kelvin
volume = R*ttotal/Ptotal;


cstar = ((gamma+1)/2)^((gamma+1)/(2*(gamma-1)))*sqrt(R*ttotal/gamma);
Astar = massflow/Ptotal*cstar;

isp = 300;
g = 9.81;
thrust = isp*massflow*g;
p_e = 101325*exp(6556.634183820265/8013.2);

syms Me
eqn = p_e/Ptotal == 1/(1+(gamma-1)/(2)*Me^2)^(gamma/(gamma-1));
Me = double(abs(vpasolve(eqn,Me)));

Tstar = ttotal*(1+((gamma-1)/2))^-1;
astar = sqrt(gamma*R*Tstar);

rhostar = massflow/(astar*Astar);

Pstar = Ptotal*(1+((gamma-1)/2))^(-gamma/(gamma-1));

rhototal = 1/((Pstar/Ptotal)^(1/gamma)/rhostar);

phi = 1;
r = m_o2/m_ch4;

Constants = ["M_p (kg)"; "Mass Flow Rate (kg/s)"; "ISP (s)"; "R (J/kg*K)"; "Gamma"; "MW (kg/kmol)"; "Equivalence Ratio"; "Mixture Ratio";"Thrust (kN)"];
Value = [mp; massflow; isp; R; gamma; mwmix; phi; r; thrust;];
ConstantValues = table(Constants, Value)



Stagnation = ["Pressure (Pa)"; "Temperature (k)"; "Chamber Volume (m3)"; "Density (kg/m3)"; ];
Value = [Ptotal; ttotal; volume; rhototal;];
stagvalues = table(Stagnation, Value)

Throat = ["Pressure (Pa)"; "Temperature (k)"; "Density (kg/m3)"; "Area (m2)"; "Mach"; "Speed of Sound (m/s)"];
Value = [Pstar; Tstar; rhostar; Astar; 1; astar;];
throatvalues = table(Throat, Value)

Texit = ttotal*(1+(gamma-1)/2*Me^2)^-1;
aexit = sqrt(gamma*R*Texit);
Aexit = Astar*((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*Me^2)^((gamma+1)/(2*(gamma-1)))/Me;
rhoexit = 1/((Pstar/p_e)^(1/gamma)/rhostar);

Exit = ["Pressure (Pa)"; "Temperature (k)"; "Density (kg/m3)"; "Area (m2)"; "Mach"; "Speed of Sound (m/s)"];
Value = [p_e; Texit; rhoexit; Aexit; Me; aexit;];
exitvalues = table(Exit, Value)






