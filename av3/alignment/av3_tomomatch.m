function av3_scan_angles_exact(refilename, PSFfilename, motlfilename, particlefilename, startindx, iterations,angincr, ...
    angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin, dispflag)
% AV3_SCANANGLES_EXACT aligns 3D subtomograms to reference
%
%   av3_scan_angles_exact(refilename, motlfilename, particlefilename, startindx, iterations,angincr,...
%       angiter,mask,hipass,lowpass,nfold,threshold,wedgelist,iclass, ibin, dispflag)
%   Filenames are expected as:
%       'filename'_#no.em
%
%   In this routine, the reference is rotated into the proposed orientation
%   (given in the motl). A X-correlation with the particle under scrutiny is
%   performed for this orientation and the ANGITER*ANGINCR neighbouring
%   orienations. 
%   A clever sampling on the unit sphere is done for psi and theta: rings around 
%   (theta_old, psi_old) are drawn. 
%   The routine takes into account the missing wedge - the semi angle of
%   the wedge SEMIANGLE has to be given. ROI is the radius if interest
%   which is used for masking the particles and RSMOOTH is the smoothing of
%   this sphere (additional to ROI). Only particles of the class ICLASS
%   (column 20 of MOTL) are correlated. For determination of CCF-peaks
%   spline interpolation is used - sub-pixel accuracy.
%
%  PARAMETERS
%   refilename          filename of reference(s) - 'refilename'_#no.em
%   PSFfilename         filename of weighting function of previous
%                       iteration
%   motlfilename        filename of corresponding motl - 'motlfilename'_#no.em
%   particlefilename    filename of 3D-particles to be aligned and averaged
%                           'particlefilename'_#no.em
%   startindx           start index - Index of first reference AND motl
%   iterations          number of iterations to be performed -after each
%                           iteration motl and average are stored
%   angincr             angular increment
%   angiter             iterations for each angle phi, psi, theta
%   mask                mask - make sure dims are all right!
%   hipass              hipass - for X-corr
%   lowpass             lowpass - for X-corr
%   nfold               symmetry of particle
%   threshold           Threshold*mean(ccc) is cutoff for averaging - e.g.
%                           choose 0.5
%   wedgelist           array containing the tilt range - 1st: no of
%                           tomogram, 2nd: minimum tilt, 3rd: max tilt
%   iclass              class of particles - default:0
%   ibin                binning - default 0
%   dispflag            flag for display - if set to 'nodisp' then no
%                       display, e.g. for batch mode - otherwise display
%
% Format of MOTL:
%    The following parameters are stored in the matrix MOTIVELIST of dimension (20, NPARTICLES):
%    column 
%       1         : Cross-Correlation Coefficient
%       2         : x-coordinate in full tomogram
%       3         : y-coordinate in full tomogram
%       4         : particle number
%       5         : running number of tomogram - used for wedgelist
%       6         : index of feature in tomogram (optional)
%       8         : x-coordinate in full tomogram
%       9         : y-coordinate in full tomogram
%       10        : z-coordinate in full tomogram
%       11        : x-shift in subvolume - AFTER rotation of template
%       12        : y-shift in subvolume - AFTER rotation of template
%       13        : z-shift in subvolume - AFTER rotation of template
%     ( 14        : x-shift in subvolume - BEFORE rotation of template )
%     ( 15        : y-shift in subvolume - BEFORE rotation of template )
%     ( 16        : z-shift in subvolume - BEFORE rotation of template )
%       17        : Phi (in deg)
%       18        : Psi
%       19        : Theta 
%       20        : class no
%
%   09/18/03 FF
%last change 19/07/04 FF - corrected bug in normalization - now CCC=1 for
%                           perfect correlation

error(nargchk(13,18,nargin))
if nargin < 15 
    iclass = 0;
end;
if nargin < 16
    ibin = 0;
end;
if nargin < 18
    dispflag = 'disp';
end;
if ibin > 0
    mask = tom_bin(mask,ibin);
end;
npixels = sum(sum(sum(mask)));
cent= [floor(size(mask,1)/2)+1 floor(size(mask,2)/2)+1 floor(size(mask,3)/2)+1];
scf = size(mask,1)*size(mask,2)*size(mask,3);
ind = startindx;


%wedge=tom_wedge(ref.Value,semiangle);

for ind = startindx:startindx+iterations-1
    name = [PSFfilename '_' num2str(ind) '.em'];
    wei_old = tom_emread(name); % GS to avoid auto correlation
    wei_old = wei_old.Value;
    name = [refilename '_' num2str(ind) '.em'];
    ref = tom_emread(name);
    disp(['read in file ' name]);
    ref = ref.Value;
    if ibin > 0
        ref = tom_bin(ref,ibin);
    end;
    average = ref*0;
    lowp = floor(size(average,1)/2)-3;%lowpass
    wei = zeros(size(average,1),size(average,2),size(average,3));%weighting function
    %mask = tom_spheremask(ones(size(average,1),size(average,2),size(average,3)),roi,rsmooth);
    name = [motlfilename '_' num2str(ind) '.em'];
    motl = tom_emread(name);
    motl = motl.Value;
    if ibin > 0
        motl(11:16) = motl(11:16)/(2^ibin);%take binning into account
    end;
    indx = find ((motl(20,:) ==1 ) | (motl(20,:) == 2) | (motl(20,:) == iclass) ); meanv = mean(motl(1,indx));
    indx = find (motl(1,:) > threshold*meanv);
    itomo_old = 0;
    for indpart = 1:size(motl,2)
        if ((motl(20,indpart) == 1) | (motl(20,indpart) == 2) | (motl(20,indpart) == iclass))
            itomo = motl(5,indpart);
            if itomo_old ~= itomo %wedge stuff - exact weighting according to list
                xx = find(wedgelist(1,:) == itomo);
                minangle= wedgelist(2,xx);
                maxangle= wedgelist(3,xx);
                wedge = av3_wedge(ref,minangle,maxangle);
                itomo_old = itomo;
            end;
            tshift = 0;
            phi_old=motl(17,indpart);
            psi_old=motl(18,indpart);
            the_old=motl(19,indpart);
            % read shift BEFORE rot
            xshift = motl(14,indpart);
            yshift = motl(15,indpart);
            zshift = motl(16,indpart);ccc = -1;
            ifile = motl(4,indpart);
            name = [particlefilename '_' num2str(ifile) '.em'];
            particle = tom_emread(name);particle = particle.Value;
            if ibin > 0
                particle = tom_bin(particle,ibin);
            end;
            particle = tom_limit(particle,-3,4,'z'); % throw away the gold
            particle4av = particle;
            tshift(1) = motl(11,indpart)/(2^ibin);
            tshift(2) = motl(12,indpart)/(2^ibin);
            tshift(3) = motl(13,indpart)/(2^ibin);
            temp_part=double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_old,-phi_old,-the_old]));
            wref = ref-real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(temp_part)).*wei_old,lowp))));
            if nfold>1
                wref = tom_symref(wref,nfold);
            end
            [mref xx1 xx2 mstd] = tom_dev(wref,'noinfo');
            wref = (wref - mref)./mstd;
            %ref = tom_limit(ref,-3,4,'z'); % throw away the gold
            % do not shift particle but mask !
            rshift(1) = motl(11,indpart);
            rshift(2) = motl(12,indpart);
            rshift(3) = motl(13,indpart);
            %rshift = [xshift yshift zshift];
            rmask = double(tom_rotate(mask,[phi_old,psi_old,the_old]));
            %rshift = tom_pointrotate(rshift,[phi_old,psi_old,the_old]);
            shiftmask = tom_shift(rmask,rshift);
            particle=shiftmask.*particle;
            particle= particle - shiftmask.*(sum(sum(sum(particle)))/npixels);%subtract mean in sphere
            fpart=fftshift(tom_fourier(particle));
            %apply bandpass
            fpart= ifftshift(tom_spheremask(fpart,lowpass,3) - tom_spheremask(fpart,hipass,2));
            %normalize
            fpart(1,1,1)=0;
            fpart = (size(fpart,1)*size(fpart,2)*size(fpart,3))*fpart/sqrt((sum(sum(sum(fpart.*conj(fpart))))));
            for phi = phi_old-angiter*angincr:angincr:phi_old+angiter*angincr
                for ithe =  0:ceil(angiter/2)
                    if ithe == 0
                        npsi=1;
                        dpsi=360;
                    else
                        %sampling for psi and the on unit sphere in rings
                        dpsi=angincr/sin(ithe*angincr/180*pi);
                        npsi = ceil(360/dpsi);
                    end;
                    for ipsi = 0:(npsi-1)
                        r = [ 0 0 1];
                        r = tom_pointrotate(r,0,ipsi*dpsi,ithe*angincr);
                        r = tom_pointrotate(r,0,psi_old,the_old);
                        the = 180/pi*atan2(sqrt(r(1).^2+r(2).^2),r(3));
                        psi = 180/pi*atan2(r(2),r(1)) + 90;
                        rref=double(tom_rotate(wref,[phi,psi,the]));
                        rref=tom_ifourier(ifftshift(tom_spheremask(wedge.*fftshift(tom_fourier(rref)))));%changed back F
                        rmask=double(tom_rotate(mask,[phi,psi,the]));
                        rref = rref.*rmask;%mask with smoothened edges
                        rref = rref - rmask.*(sum(sum(sum(rmask.*rref)))/npixels); %subtract mean in mask
                        fref=fftshift(tom_fourier(rref));
                        %apply bandpass
                        fref=ifftshift(tom_spheremask(fref,lowpass,3) - tom_spheremask(fref,hipass,2));
                        fref(1,1,1)=0;
                        %calculate rms - IMPORTANT! - missing wedge!!!
                        % changed back FF
                        % fref = fref/sqrt((sum(sum(sum(fftshift(fref.*conj(fref)).*wedge)))));
                        fref = (size(fref,1)*size(fref,2)*size(fref,3))*fref/sqrt((sum(sum(sum(fref.*conj(fref))))));% to be changed?
                        ccf = tom_spheremask(real(fftshift(tom_ifourier(fpart.*conj(fref)))),size(average,1)/5,size(average,1)/16);
                        ccf = ccf/(size(ccf,1).^3);% added for normalization - FF
                        %[pos ccctmp] = peak_det_2(real(ccf));
                        [pos ccctmp] = tom_peak(real(ccf));
                        if ccctmp > ccc
                            ccc = ccctmp;
                            ccf_max = ccf;% store CCF for best angle - idea rick gaudette
                            phi_opt=phi;
                            psi_opt=psi;
                            the_opt=the;
                            tshift = pos-cent;
                            if (strcmp(dispflag,'nodisp')~= 1 )
                                if (size(ccf)>100)
                                    binf = int8(log2(size(ccf)/64));
                                    tom_dspcub(tom_bin(ccf),binf(1));drawnow;
                                else
                                    tom_dspcub(ccf);drawnow;
                                end;
                            end;
                        end;
                    end;
                end;
            end;
            motl(17,indpart)=phi_opt;
            motl(18,indpart)=psi_opt;
            motl(19,indpart)=the_opt;
            [pos ccc] = peak_det_2(real(ccf_max));% final shift by interpolation
            tshift = pos-cent;
            motl(11,indpart) = tshift(1)*(2^ibin);
            motl(12,indpart) = tshift(2)*(2^ibin);
            motl(13,indpart) = tshift(3)*(2^ibin);
            rshift = tom_pointrotate(tshift,-psi_opt,-phi_opt,-the_opt);
            motl(14,indpart) = rshift(1)*(2^ibin);
            motl(15,indpart) = rshift(2)*(2^ibin);
            motl(16,indpart) = rshift(3)*(2^ibin);
            motl(1,indpart)=ccc;
            % take care: particle4av is NOT pre-shifted
            if (ccc > threshold*meanv) %kick off bad particles
                average = average + double(tom_rotate(tom_shift(particle4av,-tshift),[-psi_opt,-phi_opt,-the_opt]));
                if size(average,1)>100
                    tom_dspcub((tom_bin(average)));drawnow;
                else
                    tom_dspcub((average));drawnow;
                end;
                %weighting - avoid interpolation artefacts
                tmpwei = 2*tom_limit(tom_limit(double(tom_rotate(wedge,[-psi_opt,-phi_opt,-the_opt])),0.5,1,'z'),0,0.5);
                wei = wei + tmpwei;
                motl(20,indpart)=1;%good particles -> class one
            else
                motl(20,indpart)=2;%bad CCF: kick into class 2
            end;
            disp(['Particle no ' num2str(ifile) ' , Iteration no ' num2str(ind)]);
            if ( rem(indpart,10) == 0 )
                name = [motlfilename '_tmp_' num2str(ind+1) '.em'];
                tom_emwrite(name,motl);
            end;
        end; %endif 
    end;% end particle loop
    name = [motlfilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,motl);
    % do weighting
    wei = 1./wei;rind = find(wei > 100000);wei(rind) = 0;% take care for inf
    average = real(tom_ifourier(ifftshift(tom_spheremask(fftshift(tom_fourier(average)).*wei,lowp))));
    if ibin > 0
        average = tom_zoom(average,ibin);
    end;
    name = [refilename '_' num2str(ind+1) '.em'];
    tom_emwrite(name,average);
    disp(['wrote reference ' name]);
    name = [PSFfilename '_' num2str(ind+1) '.em'];    
    tom_emwrite(name,wei);
    
end; % end iteration loop