function vectsRot=tom_rotVectByAng(vect,angs,rotFlav,display,color)
%tom_rotVectByAng performs rotation of a vector for matrix of angles
%
%   vectsRot=tom_rotVectByAng(vect,ang,rotFlav,display)
%
%PARAMETERS
%
%  INPUT
%   vect                vector to be rotated (x,y,z)
%   angs                angles around the vector should be rotated
%   rotFlav            ('zxz') or 'zyz'
%   display            ('none') 'vector' or 'points' 
%   color               ([0 0 1]) colour for plotting
%
%  OUTPUT
%   binned              binned image
%
%EXAMPLE
%   
%  vectsRot=tom_rotVectByAng([1 0 0],[10 0 0; 20 0 0; 30 0 0],'zxz','vector',[0 1 0]);
%   
%
%REFERENCES
%
%SEE ALSO
%
%   created by FB 
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

if (nargin <3)
    rotFlav='zxz';
end;

if (nargin <4)
    display='none';
end;

if (nargin <5)
    color=[0 0 1];
end;

for i=1:size(angs,1)
    if (strcmp(rotFlav,'zyz'))
         angTmp=tom_eulerconvert_xmipp(angs(i,1),angs(i,2),angs(i,3));
    else
         angTmp=angs(i,:);
    end;
    vectsRot(i,:)=tom_pointrotate(vect,angTmp(1),angTmp(2),angTmp(3));
end;


if (strcmp(display,'vector'))
    ori=zeros(size(vectsRot));
    quiver3(ori(:,1),ori(:,2),ori(:,3),vectsRot(:,1),vectsRot(:,2),vectsRot(:,3),1,'color',color);
end;

if (strcmp(display,'points'))
    plot3(vectsRot(:,1),vectsRot(:,2),vectsRot(:,3),'ro','color',color);
end;


