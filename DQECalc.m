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

function DQE = DQECalc( gam, r, p, N, sensingModel, noiseModel )

% sensingModel: 1: Conventional single-pixel, 2: Scanning for single-pixel, 3: Sparse ADF-STEM, 4: ADF-STEM, each pixel has a chance p to be on, 5: Scanning ADF-STEM

P = p * N;

switch sensingModel
    case 1
        if noiseModel == 1
            DQE = ( 1 - p ) / N;
        elseif noiseModel == 2
            DQE = p .* ( 1 - p ) ./ ( P + r .* gam );
        elseif noiseModel == 3
            DQE = p .* ( 1 - p ) ./ ( gam .* r );
        end
    case 2
        if noisemodel == 1
            DQE = 1 ./ N ;
        elseif noiseModel == 2
            DQE = 1 ./ ( N .* ( 1 + gam ) );
        elseif noiseModel == 3
            DQE = 1 ./ ( N * gam );
        end
    case 4
        if noiseModel == 1
            DQE = ( 1 - p ) ./ P .* r;
        elseif noiseModel == 2
            DQE = ( 1 - p ) ./ P .* r ./ ( 1 + gam );
        elseif noiseModel == 3
            DQE = ( 1 - p ) ./ P .* r ./ gam;
        end
    case 5
        if noiseModel == 1
            DQE = r;
        elseif noiseModel == 2
            DQE = r ./ ( 1 + gam ./ r );
        elseif noiseModel == 3
            DQE = r.^2 ./ gam;
        end
    case 6
        if noiseModel == 1
            DQE = ( 2/3 - p/2 ) / N;
        elseif noiseModel == 2
            DQE = p .* ( 1/3 - p/4 ) ./ ( P/2 + r .* gam );
        elseif noiseModel == 3
            DQE = p .* ( 1/3 - p/4 ) ./ ( gam .* r );
        end
end