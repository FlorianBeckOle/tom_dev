function [avgRotEuler,avgRotMat]=tom_average_rotations(rotations,rotFlav)
%TOM_AVERAGE_ROTATIONS calculates the averages angle
%
%
%    [avgRotEuler,avgRotMat]=tom_average_rotations(rotations,rotFlav)
%   
%
%PARAMETERS
%
%  INPUT
%   rotations      nx3 matrix of rotations
%   rotFlav          ('zxz') or zyz
%
%   
%  OUTPUT
%
%   avgRotEuler      average euler angle (zxz)
%   avgRotMat         average rotation matrix
%                      
%
%EXAMPLE
%   [avgRotEuler,avgRotMat] = tom_average_rotations([10 0 0; 20 0 0]);
%
%
%REFERENCES
%
%SEE ALSO
%   tom_sum_rotation
%
%   created by FB 04/05/07
%     
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom





if (nargin<2)
    rotFlav='zxz';
end;

AngsQuat=tom_eulerconvert_Quaternion(rotations,rotFlav);
A=zeros(4,4);
M=size(AngsQuat,1);
for i=1:M
    q = AngsQuat(i,:)';
    A=(q*q')+A; 
end
[Qavg, ~]=eigs(A,1);
avgRotMat=tom_quaternion2rotMatrix(Qavg);
avgRotEuler=tom_rotmatrix2angles(avgRotMat);