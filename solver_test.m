% solver.m - Used to generate kappa coefficients to be used by kappa_analysis.py
% Ryan Gloekler, 2022
% ECSyD Laboratory, Dr. Mahdi Nikdast
% Last update April 8, 2022

w = 800;
R = 10000;
d = 100;
kappa = 0.0567189;

syms ae ao ge go real
eqn = (pi / 1550) * (ae/ge * (exp(-1 * ge * d) * sqrt(2 * pi * ge*(R + w/2))) + ao/go*exp(-1 * go * d) * sqrt(2*pi*go*(R + w/2)))== kappa;


varlimits = [0 0.3; 0 0.3; 0 0.3; 0 0.3];
[s1, s2, s3, s4] = vpasolve(eqn, ae, ao, ge, go, varlimits);


answer = (pi / 1550)*(s1/s3*exp(-s3*100)*sqrt(2*pi*s3*(10000 + 400/2)) + s2/s4*exp(-s4*100)*sqrt(2*pi*s4*(10000 + 400/2)));

disp(s1);
disp(s2);
disp(s3);
disp(s4);
