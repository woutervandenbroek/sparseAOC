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

function Is = fisherInformation( s, As, noiseModel, I0, c )

switch noiseModel
    case 1  % Pure Poisson noise
        [ ~, n ] = size( As );
        y = As * s;
        Atmp = As ./ repmat( y, [ 1, n ] );
        Is = I0 * ( As' * Atmp );
        clear Atmp
        
    case 2 % Poisson + read-out noise
        [ ~, n ] = size( As );
        y = As * s;
        Atmp = As ./ repmat( y + c / I0, [ 1, n ] );
        Is = I0 * ( As' * Atmp );
        Atmp = 0.5 * ( Atmp ./ repmat( y + c / I0, [ 1, n ] ) );
        Is = Is + As' * Atmp;
        clear Atmp
        
    case 3 % Only read-out noise
        Is = I0^2 * ( As' * As ) / c;
end