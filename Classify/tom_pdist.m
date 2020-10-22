function [dists]=tom_pdist(in,dmetric,in_Inv)
%TOM_PDIST calculates distances between observations analog to pdist2
%
%    dists=tom_pdist(in,dmetric,in_Inv)
%
%PARAMETERS
%
%  INPUT
%   in                 inputdata nxdimension for angles nx3 in zxz 
%   dmetric        ('euc') for euclidean distance metric or  
%                        'ang'   for angles
%   in_Inv           ('')  inverse data neede for needed for transformations    
%   verbose        (1) 0 for no output
%
%  OUTPUT
%   dists             distances in the same order as pdist from matlab
%
%EXAMPLE
%  
% parpool('local',32); 
% dd=tom_pdist([0 0 0;0 0 10; 10 20 30],'ang');
% figure; imagesc(dd); colormap gray;
%
%REFERENCES
%
%SEE ALSO
%   pdist
%
%   created by FB 03/14/19
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
    dmetric='euc';
end;

if (nargin<3)
    in_Inv='';
end;


verbose=1;

jobList=genJobList(size(in,1));

%dists=zeros(size(jobList,1),1,'single');
%  inOrg=in;
%  disp('alloc')
%  in=gpuArray(inOrg);
% disp('alloc done');


if (strcmp(dmetric,'euc'))
    g1=in(jobList(:,1),:);
    g2=in(jobList(:,2),:);
    dv=g2-g1;
    dists=sqrt(sum(dv.*dv,2));
    if (isempty(in_Inv)==0)
          g1Inv=in_Inv(jobList(:,1),:);
          dv=g2-g1Inv;
          distsInv=sqrt(sum(dv.*dv,2));
          dtmp=min([dists distsInv],[],2);
          dists=dtmp;
    end;
end;

if (strcmp(dmetric,'ang'))
    [dists,RinInv]=calcAngDist(in,jobList);
     if (isempty(in_Inv)==0)
        distsInv=calcAngDist(in_Inv,jobList,RinInv);
        dtmp=min([dists distsInv],[],2);
        dists=dtmp;   
     end;
end;

dists=dists';



function [dists,Rsback]=calcAngDist(in,jobList,Rcmb)

Rsback='';
if (nargin<3)
   Rcmb='';
end;


Rin=zeros(3,3,size(in,1));
 if (isempty(Rcmb))
       RinInv=zeros(3,3,size(in,1));
end;



wb=tom_progress(size(in,1),'calc rot matrices');
%parfor i=1:size(in,1) %does not work for invers no reason why
for i=1:size(in,1)
    [~,~,Rin(:,:,i)]=tom_sum_rotation(in(i,:),[0 0 0]);
    if (isempty(Rcmb))
       RinInv(:,:,i)=inv(Rin(:,:,i));
    end;
    if (mod(i,100))
        wb.update();
    end;
end;
clear('wb');
Rs=Rin(:,:,jobList(:,1));
disp('calc matrix done');

if (isempty(Rcmb))
    RsInv=RinInv(:,:,jobList(:,2));
else
    RsInv=Rcmb;
end;


Rp=tom_multiprod(Rs,RsInv);
Reye=repmat(eye(3,3),1,1,size(Rs,3));
Rp=Rp.*Reye;
Rp=tom_multiprod(Rp,[1 1 1]');
dists=sum(Rp,1);
dists=squeeze(dists);
dists=acosd((dists-1)./2);

if (exist('RinInv','var'))
    if (isempty(RinInv)==0)
        Rsback=RinInv(:,:,jobList(:,2));
    end;
end;

function jobList=genJobList(szIn)

for i=1:szIn
    v2=(i+1):szIn;
    v1=ones(length(v2),1)*i;
    jobList{i}=[v1 v2'];
end;
jobList=cell2mat(jobList');


