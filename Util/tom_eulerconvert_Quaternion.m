function Q=tom_eulerconvert_Quaternion(angels,rotFlav)
%TOM_EULERCONVERT_QUATERNION converts a euler angles to quaternions
%
%   Q=tom_eulerconvert_Quaternion(angels,rotFlav)
%
%   tom_eulerconvert_Quaternion converts a matrix of euler ang to a matrix
%   of Quaternions
%
%PARAMETERS
%
%  INPUT
%   angels                 nx3 matrix of euler angles
%   rotFlav                 ('zxz') or 'zyz'
%
%  OUTPUT
%   euler_out            resulting Euler angles
%
%EXAMPLE
%   [euler_out]=tom_eulerconvert_Quaternion([10 20 30; 0 0 20]);
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 24/08/19
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

if (nargin< 2)
    rotFlav='zxz';
end

Q=zeros(size(angels,1),4);

for ii=1:size(angels,1)
    
    if  (strcmp(rotFlav,'zyz'))
       [~,angTmp]=tom_eulerconvert_xmipp(angels(ii,:));
    else
        angTmp=angels(ii,:);
    end;
    
    [~,~,R]=tom_sum_rotation(angTmp,[0 0 0]);
    
    [r,c] = size( R );
    if( r ~= 3 || c ~= 3 )
        fprintf( 'R must be a 3x3 matrix\n\r' );
        return;
    end;
    
    Rxx = R(1,1); Rxy = R(1,2); Rxz = R(1,3);
    Ryx = R(2,1); Ryy = R(2,2); Ryz = R(2,3);
    Rzx = R(3,1); Rzy = R(3,2); Rzz = R(3,3);
    
    w = sqrt( trace( R ) + 1 ) / 2;
    
    % check if w is real. Otherwise, zero it.
    if( imag( w ) > 0 )
        w = 0;
    end;
    
    x = sqrt( 1 + Rxx - Ryy - Rzz ) / 2;
    y = sqrt( 1 + Ryy - Rxx - Rzz ) / 2;
    z = sqrt( 1 + Rzz - Ryy - Rxx ) / 2;
    
    [~, i ] = max( [w,x,y,z] );
    
    if( i == 1 )
        x = ( Rzy - Ryz ) / (4*w);
        y = ( Rxz - Rzx ) / (4*w);
        z = ( Ryx - Rxy ) / (4*w);
    end;
    
    if( i == 2 )
        w = ( Rzy - Ryz ) / (4*x);
        y = ( Rxy + Ryx ) / (4*x);
        z = ( Rzx + Rxz ) / (4*x);
    end;
    
    if( i == 3 )
        w = ( Rxz - Rzx ) / (4*y);
        x = ( Rxy + Ryx ) / (4*y);
        z = ( Ryz + Rzy ) / (4*y);
    end;
    
    if( i == 4 )
        w = ( Ryx - Rxy ) / (4*z);
        x = ( Rzx + Rxz ) / (4*z);
        y = ( Ryz + Rzy ) / (4*z);
    end;
    
    Q(ii,:) = [ w; x; y; z ];
    
end;






