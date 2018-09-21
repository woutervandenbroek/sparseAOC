% Copyright 2018 Wouter Van den Broek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% sparseAOC lets you calculate A-optimality criteria (AOCs) for mcOuter 
% different settings. Since the AOC is a statistical variable itself, one 
% can repeat the calculation mcInner times and get the average and the 
% standard deviation.

% This is the software used to calculate the data for Figures 3 through 6
% in the paper [REFERENCE]

clear
close all

% Number of iterations
% Number of Monte Carlo realizations of each of the mcOuter settings
mcInner = 1;
% Number of different settings
mcOuter = 10;

% 1 for piecewise linear and 0 for piecewise constant test object
flagLin = 1; 

% Size of the n x n test object
n = 100;
N = n * n;
% What sensing matrix to use
% sensingModel 1: Single-pixel with 0/1 pixels
%              2: Single-pixel with fractional pixels
%              3: ADF-STEM CS, p on-pixels
%              4: ADF-STEM, scanning
sensingModel = 1;
% noiseModel 1: Pure Poisson noise
%            2: Poisson + read-out noise
%            3: only read out noise
noiseModel = 2;
% Electron or photon dose y = I0 A x
I0 = 7.460e5;

checkConsistency( sensingModel, noiseModel );

% Set the upper and lower bound of the relative read-out noise
gamLo = 0.99999 * 625;
gamHi = 1.00001 * 625;
% Set the upper and lower bound of the fraction of on-pixels
pLo = 0.01;
pHi = 0.50;

% Make the phantom
if flagLin == 0  % For piecewise constant test object
    x = phantom( n );
    x( x < 0 ) = 0;
    % Add small background to mimic real experiments better
    x = x * 0.9 + 0.1; % 0.1 < x < 1
end
if flagLin == 1  % For piecewise linear test object
    load('rampDiscs100x100.mat')
    n = 100;
    N = n * n;
end

figure(1);
subplot(2,2,1); imagesc( x )  ; colormap gray, axis equal, axis tight, colorbar
title( 'Test object' )

% Find the sparse "basis" for the 1st direction and the transformation 
% matrices for x = DinvS * s and t = DT * x
if flagLin == 0
    [s, t, DT, DinvS ] = sparseComponents( x );
end
if flagLin == 1
    [ s, DinvS ] = sparseComponentsLinear( x );
end
K = numel( s );

% Everything in long-vector format
x = x(:);

% Set the upper and lower bound of the reduction factor; M = round( r * N )
rLo = 2000.1 / N;
rHi = 1999.9 / N;

% (Randomly) initialize gam, r and p
gamAr = rand( mcOuter, 1 ) * ( gamHi - gamLo ) + gamLo ;
rAr   = rand( mcOuter, 1 ) * ( rHi - rLo ) + rLo;
% Array of fraction of on-pixels:
% pAr   = linspace( pLo, pHi, mcOuter )';  % Linearly spaced
pAr = exp( linspace( log( pLo ), log( pHi ), mcOuter )' );  % Logarithmically spaced

% The AOCs of s (ST) and x (X).  mu is the average AOC, std the standard
% deviation on the AOC.  Due to legacy, AOC is denoted CRB in the 
% following.
CRBSTmu = zeros( mcOuter, 1 );  
CRBSTstd = zeros( mcOuter, 1 );
CRBXmu = zeros( mcOuter, 1 );
CRBXstd = zeros( mcOuter, 1 );

reverseStr = '';
cnt = 0;

for j = 1:mcOuter
    % Primitive progress counter
    msg = sprintf('   Progress: %i pct.', round( cnt / mcOuter * 100 ) );
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    cnt = cnt + 1;
    
    % Copy gam, r and p out of the array
    gam = gamAr( j );
    r   = rAr( j );
    p   = pAr( j );

    % Translate those to the other variables
    M = round( r * N );
    r = M / N;
    % Calculate the read-out noise variance
    [c, gam] = gam2C( x, I0, gam, sensingModel, noiseModel );
    % fraction of on-pixels consistent with integer value for N
    p = round( p * N ) / N;
    
    % Store for later
    gamAr( j ) = gam;
    rAr( j ) = r;
    pAr( j ) = p;

    varTV = zeros( mcInner, 1 );
    varX  = zeros( mcInner, 1 );
    % Actual simulations start
    for k = 1:mcInner        
        As = sensingMatrix( x, M, sensingModel, p );
        As = sparse( As );
        As = As * DinvS;

        % Fisher information for s
        Is = fisherInformation( s, As, noiseModel, I0, c );

        % AOC of s
        covS = info2AOC( Is );

        % AOC of the total (generalized) variation
        if flagLin == 0
            [ varTV( k ), varX( k ) ] = meanSquaredError( covS, DT, DinvS, K );
        end
        if flagLin == 1
            varTV( k ) = trace( covS ) / K;
            varX( k ) = trace( DinvS * ( covS * ( DinvS' ) ) ) / size( DinvS, 1 );
        end
    end
    
    CRBSTmu( j )  = mean( varTV );
    CRBSTstd( j ) = std( varTV, 0 );
    CRBXmu( j )   = mean( varX );
    CRBXstd( j )  = std( varX, 0 );
end
disp( ' ' )

% Calculate the corresponding DQEs
DQE = DQECalc( gamAr, rAr, pAr, N, sensingModel, noiseModel );

% Scaling 1/DQE to CRBmu with a least absolute difference fit
alfaST = sum( CRBSTmu ./ DQE ) / sum( DQE.^(-2) );
for j = 1:100
    alfaST = ( sum( CRBSTmu ./ DQE ./ sqrt( ( CRBSTmu - alfaST * DQE.^-1 ).^2 + 1e-12 ) ) ) ./ ( sum( DQE.^-2 ./ sqrt( ( CRBSTmu - alfaST * DQE.^-1 ).^2 + 1e-12 ) ) );
end

% Produce nice figures of the result
nP = ceil( ( pHi - pLo ) / 0.002 );
pPlot = linspace( pLo, pHi, nP );
DQEPlot = DQECalc( gamAr(1), rAr(1), pPlot, N, sensingModel, noiseModel );

figure(2);
plot( pPlot, alfaST./DQEPlot, 'r', pAr, CRBSTmu, 'bo' )
xlabel( '1 / DQE' )
ylabel( 'CRBSTmu' )

% Scaling 1/DQE to CRBmu with a least absolute difference fit
alfaX = sum( CRBXmu ./ DQE ) / sum( DQE.^(-2) );
for j = 1:100
    alfaX = ( sum( CRBXmu ./ DQE ./ sqrt( ( CRBXmu - alfaX * DQE.^-1 ).^2 + 1e-12 ) ) ) ./ ( sum( DQE.^-2 ./ sqrt( ( CRBXmu - alfaX * DQE.^-1 ).^2 + 1e-12 ) ) );
end

figure(3);
plot( pPlot, alfaX./DQEPlot, 'r', pAr, CRBXmu, 'bo' )
xlabel( '1 / DQE' )
ylabel( 'CRBXmu' )
