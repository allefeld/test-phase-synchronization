function errorsFig

% Code for Allefeld & Kurths, "Testing for Phase Synchronization"
% checking validity of tests in a simulation
% generates figure using results of testing (testing.mat)
%
% Copyright (C) 2003 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


    load('testing')
    
    subplot(3, 1, 1)
    sa = sqrt(0.05 * (1 - 0.05) / N);
    plot(rhos', sH0', 'x-', [0 1], [0.05 0.05-sa 0.05+sa]' * [1 1], 'k-')
    axis([0 1 0 0.08])
    legend({'Parametric', 't Direct Variance', 't Bootstrap Variance', 'Bootstrap H0', 'Permutation H0'}, 4)
    set(gca, 'FontSize', 12)
    xlabel('\rho', 'FontSize', 14)
    ylabel('Error of the first kind', 'FontSize', 14)
    title('a)', 'FontSize', 18)
    
    subplot(3, 1, 2)
    plot(rhos', sH1a', '-')
    set(gca, 'FontSize', 12)
    xlabel('\rho_1', 'FontSize', 14)
    ylabel('Power function', 'FontSize', 14)
    title('b)', 'FontSize', 18)
    
    subplot(3, 1, 3)
    plot(rhos', sH1b', '-')
    hold on
    plot(rhos', sH1c', '-')
    set(gca, 'FontSize', 12)
    xlabel('\rho_1', 'FontSize', 14)
    ylabel('Power function', 'FontSize', 14)
    title('c)', 'FontSize', 18)

    