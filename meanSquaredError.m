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


function [ mseST, mseX ] = meanSquaredError( covS, DT, DinvS, K )

mseST = trace( covS );
DT = DT * DinvS;
mseST = mseST + trace( DT * ( covS * ( DT' ) ) );
mseST = 0.5 * mseST / K;

mseX = trace( DinvS * ( covS * ( DinvS' ) ) ) / size( DinvS, 1 );