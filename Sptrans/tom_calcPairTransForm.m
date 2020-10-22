function [posTr,angTr,lenPosTr,lenAngTr]=tom_calcPairTransForm(pos1,ang1,pos2,ang2,dMetric)
%TOM_CALCPAIRTRANSFORM calculates rel transformaton between two poses
%             
%
%    [posTr,angTr]=tom_calcPairTransForm(pos1,ang1,pos2,ang2)
%
%PARAMETERS
%
%  INPUT
%   pos1                 position 1
%   ang1                 angle1 in zxz
%   pos2                 posion2
%   ang2                 angle2 in zxz
%   dMetric             ('exact') at the moment only exact implemented
%
%  OUTPUT
%   posTr              transformation vector between two points 
%   angTr              transformation angle between two points
%   lenPosTr         length of transformation vector
%   lenAngTr        angular distance from [0 0 0] to angTr
%
%
%EXAMPLE
%   [pos,rot]=tom_calcPairTransForm([1 1 1],[0 0 10],[2 2 2],[0 0 30]);
%  
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 12/08/19
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
%

if (strcmp(dMetric,'exact'));
    
    %rotate 2 by the inverse of 1
    ang1Inv=[-ang1(2) -ang1(1) -ang1(3)];
    angTr=tom_sum_rotation([ang2;ang1Inv],[0 0 0;0 0 0]);
    
    pos2Rel=pos2-pos1;
    posTr=tom_pointrotate(pos2Rel,ang1Inv(1),ang1Inv(2),ang1Inv(3));
    
    if (nargout>2)
        lenPosTr=norm(posTr);
    end;
    
    if (nargout>3)
        lenAngTr=tom_angular_distance([0 0 0],angTr);
    end;
    
end;

