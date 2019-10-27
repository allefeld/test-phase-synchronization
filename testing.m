function testing

% Code for Allefeld & Kurths, "Testing for Phase Synchronization"
% checking validity of tests in a simulation
% results are saved to file testing.mat
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


    rhos = [eps, 0.02 :0.02: 0.98];
    N = 4000;

    fprintf('H0:')
    sH0 = [];
    for rho1 = rhos
        fprintf(' %3.3f', rho1)
        rho2 = rho1;
        sH0 = [sH0 rejprob(rho1, rho2, N)];
    end
    fprintf('\n')
    save('testing', 'rhos', 'N', 'sH0')
    
    fprintf('H1a:')
    sH1a = [];
    rho1 = eps;
    for rho2 = rhos
        fprintf(' %3.3f', rho2)
        sH1a = [sH1a rejprob(rho1, rho2, N)];
    end
    fprintf('\n')
    save('testing', 'rhos', 'N', 'sH0' ,'sH1a')

    fprintf('H1b:')
    sH1b = [];
    rho1 = 0.4;
    for rho2 = rhos
        fprintf(' %3.3f', rho2)
        sH1b = [sH1b rejprob(rho1, rho2, N)];
    end
    fprintf('\n')
    save('testing', 'rhos', 'N', 'sH0' ,'sH1a', 'sH1b')

    fprintf('H1c:')
    sH1c = [];
    rho1 = 0.8;
    for rho2 = rhos
        fprintf(' %3.3f', rho2)
        sH1c = [sH1c rejprob(rho1, rho2, N)];
    end
    fprintf('\n')
    save('testing', 'rhos', 'N', 'sH0' ,'sH1a', 'sH1b' ,'sH1c')


function p = rejprob(rho1, rho2, N)
% calculates rejection probabilities in five different tests,
% for a wrapped normal distribution, depending on given rhos

    n = 100;
    sigma1 = sqrt(-2 * log(rho1));
    sigma2 = sqrt(-2 * log(rho2));

    p = 0;
    for i = 1 : N

        sig(1) = testP(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
        sig(2) = testTDV(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
        sig(3) = testTBV(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
        sig(4) = testBH0(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
        sig(5) = testPH0(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
% alternative implementations ("parametric bootstrap"), run faster
%        sig(3) = testTBVp(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
%        sig(4) = testBH0p(sigma1 * randn(n, 1), sigma2 * randn(n, 1));
%        sig(5) = testPH0p(sigma1 * randn(n, 1), sigma2 * randn(n, 1));

        p = p + sig;
    end
    p = p(:) / N;


function sig = testP(theta1, theta2)
% parametric test for wrapped normal, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    sig = (abs(atanh(R1) - atanh(R2)) * sqrt(n) > 1.9600);
    
    
function sig = testTDV(theta1, theta2)
% non-parametric t-test based on direct variance estimation, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    C1 = mean(cos(theta1));
    S1 = mean(sin(theta1));
    R1 = sqrt(C1 ^ 2 + S1 ^ 2);
    t1 = atan2(S1, C1);
    V1 = mean((cos(theta1 - t1) - R1) .^ 2) / (n - 1);
    C2 = mean(cos(theta2));
    S2 = mean(sin(theta2));
    R2 = sqrt(C2 ^ 2 + S2 ^ 2);
    t2 = atan2(S2, C2);
    V2 = mean((cos(theta2 - t2) - R2) .^ 2) / (n - 1);
    sig = (2 * (1 - cdfStudent(2 * (n - 1), abs(R1 - R2) / sqrt(V1 + V2))) < 0.05);
    
    
function sig = testTBV(theta1, theta2)
% non-parametric t-test based on bootstrap variance estimation, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    V1 = var(bootReps(theta1, n, 200));
    V2 = var(bootReps(theta2, n, 200));
    sig = (2 * (1 - cdfStudent(2 * (n - 1), abs(R1 - R2) / sqrt(V1 + V2))) < 0.05);

    
function sig = testTBVp(theta1, theta2)
% non-parametric t-test based on parametric bootstrap variance estimation, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    V1 = var(pBootReps(theta1, n, 200));
    V2 = var(pBootReps(theta2, n, 200));
    sig = (2 * (1 - cdfStudent(2 * (n - 1), abs(R1 - R2) / sqrt(V1 + V2))) < 0.05);

    
function sig = testBH0(theta1, theta2)
% non-parametric test based on bootstrap H0 distribution simulation, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    theta0 = [theta1; theta2];
    Rdiffs = sort(abs(bootReps(theta0, n, 4000) - bootReps(theta0, n, 4000)));
    RdiffUp = mean(Rdiffs(end - 200 : end - 200 + 1));
    sig = (abs(R1 - R2) > RdiffUp);
    
    
function sig = testBH0p(theta1, theta2)
% non-parametric test based on parametric bootstrap H0 distribution simulation, alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    theta0 = [theta1; theta2];
    Rdiffs = sort(abs(pBootReps(theta0, n, 4000) - pBootReps(theta0, n, 4000)));
    RdiffUp = mean(Rdiffs(end - 200 : end - 200 + 1));
    sig = (abs(R1 - R2) > RdiffUp);
    
    
function sig = testPH0(theta1, theta2)
% permutation test (H0 distribution simulation), alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    theta0 = [theta1; theta2];
    Rdiffs = sort(abs(permReps(theta0, n, n, 4000)));
    RdiffUp = mean(Rdiffs(end - 200 : end - 200 + 1));
    sig = (abs(R1 - R2) > RdiffUp);


function sig = testPH0p(theta1, theta2)
% permutation test (H0 distribution simulation), alpha = 0.05
% assumption: theta1, theta2 column vectors of equal size
    n = size(theta1, 1);
    R1 = abs(mean(exp(i * theta1)));
    R2 = abs(mean(exp(i * theta2)));
    theta0 = [theta1; theta2];
    Rdiffs = sort(abs(pPermReps(theta0, n, 4000)));
    RdiffUp = mean(Rdiffs(end - 200 : end - 200 + 1));
    sig = (abs(R1 - R2) > RdiffUp);


function Rs = bootReps(theta, n, M)
% generates M bootstrap replications of the statistic R
% on a sample of size n from the distribution of theta
% assumption: theta column vector
    boot = floor(rand(n, M) * size(theta, 1)) + 1;
    eitheta = exp(i * theta);
    Rs = abs(mean(eitheta(boot)));
    

function Rs = pBootReps(theta, n, M)
% generates M bootstrap replications of the statistic R
% on a sample of size n from the distribution of theta
% parametrically by the binormal assumption for the distribution of (C, S)
% assumption: theta column vector
    % trigonometric moments of the distribution of theta
    a = mean(cos(theta));
    b = mean(sin(theta));
    a2 = mean(cos(2 * theta));
    b2 = mean(sin(2 * theta));
    % moments of the distribution of (C, S)
    m = [a ; b];
    C = 1 / (2 * n) * [1 + a2 - 2 * a ^ 2, b2 - 2 * a * b;
                       b2 - 2 * a * b,     1 - a2 - 2 * b ^ 2];
    % generate samples of (C, S) by linear transformation of a standard binormal
    z = randn(2, M);
%     [V, D] = eig(C);
%     d = V * sqrt(D) * z + m * ones(1, M);
    R = chol(C);
    d = R' * z + m * ones(1, M);
    % calculate values of R
    Rs = sqrt(d(1, :) .^ 2 + d(2, :) .^ 2);

    
function Rdiffs = pPermReps(theta, n, M)
    % trigonometric moments of the distribution of theta
    a = mean(cos(theta));
    b = mean(sin(theta));
    a2 = mean(cos(2 * theta));
    b2 = mean(sin(2 * theta));
    % moments of the distribution of (C, S)
    m = [a ; b];
    C = 1 / (2 * n) * [1 + a2 - 2 * a ^ 2, b2 - 2 * a * b;
                       b2 - 2 * a * b,     1 - a2 - 2 * b ^ 2];
    % generate samples of (C, S) by linear transformation of a standard binormal
    z = randn(2, M);
%     [V, D] = eig(C / 2);
%     d1 = m * ones(1, M) + V * sqrt(D) * z;
%     d2 = m * ones(1, M) - V * sqrt(D) * z;
    R = chol(C / 2);
    d1 = m * ones(1, M) + R' * z;
    d2 = m * ones(1, M) - R' * z;
    % calculate values of R
    Rdiffs = sqrt(d1(1, :) .^ 2 + d1(2, :) .^ 2) - sqrt(d2(1, :) .^ 2 + d2(2, :) .^ 2);


function Rdiffs = permReps(theta, n, m, M)
% generates M permutation replications of the statistic R1 - R2
% on samples of sizes n and m from the distribution of theta
% assumption: theta of size n + m
    [dummy, perm] = sort(rand(n + m, M));
    eitheta = exp(i * theta);
    eitheta = eitheta(perm);
    Rdiffs = abs(mean(eitheta(1 : n, :))) - abs(mean(eitheta(n + 1 : end, :)));
    

function P = cdfStudent(f, x)
% cumulative distribution function of the Student t distribution with f degrees of freedom
    P = 1/2 * (1 + (1 - betainc(f ./ (f + x .^ 2), f / 2, 1/2)) .* sign(x));
