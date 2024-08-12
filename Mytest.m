% Calculate the impedance of the dipole created in example 1 over the
% frequency span 50MHz - 100MHz.

clc;clear;close all
p = patchMicrostrip('Length',75e-3, 'Width', 37.5e-3,               ...
        'GroundPlaneLength', 120e-3, 'GroundPlaneWidth', 120e-3,        ...
        'FeedOffset', [-18.75e-3 0]);
figure()
impedance(p,linspace(1.70e9,1.75e9,51))
figure()
pattern(p, 1.75e9)
figure()
current(p, 1.75e9, "metal")
figure()
beamwidth(p, 1.75e9, 1:1:360, 140)
show(p)