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

function [ s, DinvS ] = sparseComponentsLinear( x )

n = size( x, 1 );

% Forward filter: s = D * x
D = zeros( n, n );
D( 1:2, 1:2 ) = [ -1 0.25; 0.25 0 ];
D( n, 1:2 )   = [ 0.25 0 ];
D( 1:2, n )   = [ 0.25; 0 ];

s = real( ifft2( fft2( x ) .* fft2( D ) ) );

tmp = zeros( n, n );
tmp( 2:(n-1), 2:(n-1) ) = s( 2:(n-1), 2:(n-1) );
s = tmp;

N = n * n;

s( abs( s ) < 10 * eps ) = 0;
s = s(:);

% imagesc( reshape( s, n, n ) )  ; colormap gray, axis equal, axis tight, colorbar
% title( 'Sparse components s' )

% Linear index of the non-zeros in s
vInd = find( s );
vInd = vInd(:);

% The anti-derivatives
tmp = fft2( D );
G = 1./ tmp;
G(1) = 0;
g = real( ifft2( G ) );
clear tmp G

DinvS = zeros( N, N );
for j2 = 1:n
    tmp = circshift( g, j2 - 1, 2 );
    tmp = circshift( tmp, -1, 1 );
    for j1 = 1:n
        tmp = circshift( tmp, 1, 1 );
        DinvS( ( j2 - 1 ) * n + j1, : ) = tmp(:);
    end
end

% Add the DC components and update DT so it produces the DCs too
s = [ s; mean( x(:) ) ];

DinvS = [DinvS, ones( N, 1 )];

% Remove the zeros from 
vInd = [ vInd; N+1 ];
s = s( vInd );
DinvS = DinvS( :, vInd );
























