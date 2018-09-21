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

function [s, t, DT, DinvS ] = sparseComponents( x )

n = size( x, 1 );
N = n * n;

% Horizontal derivative
s = x - x( :, [2:n, 1] );
s = s(:);
subplot(2,2,2); imagesc( reshape( s, n, n ) )  ; colormap gray, axis equal, axis tight, colorbar
title( 'Horizontal derivative' )

% Vertical derivative
t = x - x( [2:n, 1], : );
t = t(:);
subplot(2,2,3); imagesc( reshape( t, n, n ) )  ; colormap gray, axis equal, axis tight, colorbar
title( 'Vertical derivative' )

% Vertical derivative operator. Horizontal one is not needed
DT = zeros( N, N );
g = [ 1; -1; zeros( n - 2, 1 ) ];
for j2 = 1:n
    tmp = zeros( n, n );
    for j1 = 1:n
        tmp( :, j2 ) = circshift( g, [j1-1 0] );
        DT( ( j2 - 1 ) * n + j1, : ) = tmp(:);
    end    
end

% Make the variation
v = sqrt( s.^2 + t.^2 );
subplot(2,2,4); imagesc( reshape( v, n, n ) )  ; colormap gray, axis equal, axis tight, colorbar
title( 'Local variation' )

% Linear index of the non-zeros in the variation
vInd = find( v );
vInd = vInd(:);

% The anti-derivatives
h = zeros( n, 1 );
h( [1 2] ) = [1 -1];
H = fft( h );
G = 1./ H;
G(1) = 0;
g = real( ifft( G ) );
clear h H G

% Horizontal anti-derivative
DinvS = zeros( N, N );
for j1 = 1:n
    tmp = zeros( n, n );
    for j2 = 1:n
        tmp( j1, : ) = circshift( g, [j2-1 0] );
        DinvS( ( j2 - 1 ) * n + j1, : ) = tmp(:);
    end
end

% Add the DC components and update DT so it produces the DCs too
s = [ s; mean( x, 2 ) ];
t = [ t; mean( x, 1 )'];

DinvS = [DinvS, zeros( N, n )];
DT    = [DT;    zeros( n, N )];
for j = 1:n
    DinvS( ( (1:n:N) - 1 ) + j, N + j) = 1;
    DT( N + j, (1:n) + (j-1)*n ) = 1/n;
end

% Remove the zeros from v
vInd = [ vInd; (( N + 1 ):( N + n ))' ];
s = s( vInd );
DinvS = DinvS( :, vInd );
t = t( vInd );
DT = DT( vInd, : );
























