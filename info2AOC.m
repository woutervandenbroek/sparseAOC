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

function covS = info2AOC( Is )

nS = size( Is, 1 );
rnk = rank( Is );
if rnk == nS
    covS = inv( Is );
else
    disp( ['    !!! WARNING in info2CRLB: Is is rank deficient: ', int2str( rnk ), ' / ', int2str( nS ), ' !!! '] )
    [ U, S, V ] = svd( Is );
    tmp = diag( S );
    Srnk = diag( [ tmp(1:rnk).^(-1); zeros( nS - rnk,1 ) ] );
    covS = V * Srnk * U';
end