function R = tom_quaternion2rotMatrix(Qrotation)
%TOM_QUATERNION2ROTMATRIX calculates rotation matrix from Quat
%   R = tom_quaternion2rotMatrix(Qrotation)
%
%PARAMETERS
%
%  INPUT
%   Qrotation     input Quaternion
%  
%  OUTPUT
%    R                 rotation Matrix
%                          
%EXAMPLE
%   
%   R = tom_quaternion2rotMatrix([0.9330    0.2578    0.0226    0.2500]);
%  
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 29/07/19
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

w = Qrotation( 1 );
x = Qrotation( 2 );
y = Qrotation( 3 );
z = Qrotation( 4 );
Rxx = 1 - 2*(y^2 + z^2);
Rxy = 2*(x*y - z*w);
Rxz = 2*(x*z + y*w);
Ryx = 2*(x*y + z*w);
Ryy = 1 - 2*(x^2 + z^2);
Ryz = 2*(y*z - x*w );
Rzx = 2*(x*z - y*w );
Rzy = 2*(y*z + x*w );
Rzz = 1 - 2 *(x^2 + y^2);
R = [ 
    Rxx,    Rxy,    Rxz;
    Ryx,    Ryy,    Ryz;
    Rzx,    Rzy,    Rzz];

