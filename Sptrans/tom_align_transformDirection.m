function pairListAlg=tom_align_transformDirection(pairList,outputName)
%TOM_ALIGN_TRANSFORMDIRECTION aligns vector direction for given pair.star
%                                                           each transformation class is treated seperatly  
%
%   pairListAlg=tom_align_transformDirection(pairList)
%
%PARAMETERS
%
%  INPUT
%   pairList               pari star file                  
%   outputName       (opt.) name of the output pair star file
%  
%  OUTPUT
%   pairListAlg          aligned pair list
%
%EXAMPLE
%
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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


%parse inputs
if (ischar(pairList))
    pairList=tom_starread(pairList);
end;

allClasses=[pairList(:).pairClass];
allClassesU=unique(allClasses);

pairListAlg=pairList;
for i=1:length(allClassesU)
    idx=find(allClasses==allClassesU(i));
    pairListAlg(idx)=alignDir(pairList(idx));
end;

if (nargin>1)
    tom_starwrite(outputName,pairListAlg);
end;


function pairList=alignDir(pairList)


vects=[pairList(:).pairTransVectX ; pairList(:).pairTransVectY ; pairList(:).pairTransVectZ];
vects=vects';
vectsInv=[pairList(:).pairInvTransVectX ; pairList(:).pairInvTransVectY ; pairList(:).pairInvTransVectZ];
vectsInv=vectsInv';

if (size(vects,1)>1)
    cl=kmeans(vects,2);
else
    cl=1;
end;

if (length(find(ismember(cl,1)))>length(find(ismember(cl,2))) )
    useCl=1;
else
    useCl=2;
end;
idx=find(cl==useCl);
meanV=mean(vects(idx,1:3),1);


diffV=vects-repmat(meanV,size(vects,1),1);
diffV=sqrt(sum(diffV(:,:).*diffV(:,:),2));
diffVInv=vectsInv-repmat(meanV,size(vectsInv,1),1);
diffVInv=sqrt(sum(diffVInv(:,:).*diffVInv(:,:),2));

idxSwap=find(diffV<diffVInv);

for i=1:length(idxSwap)
    pairList(idxSwap(i))=swapPairOrderEntry(pairList(idxSwap(i)));
end;

vects(idxSwap,:)=vectsInv(idxSwap,:);   
meanV=mean(vects);



function entry=swapPairOrderEntry(entry)


entry=swapFields(entry,'pairIDX1','pairIDX2');

%TransVect
entry=swapFields(entry,'pairTransVectX','pairInvTransVectX');
entry=swapFields(entry,'pairTransVectY','pairInvTransVectY');
entry=swapFields(entry,'pairTransVectZ','pairInvTransVectZ');

%TransAng
entry=swapFields(entry,'pairTransAngleZXZPhi','pairInvTransAngleZXZPhi');
entry=swapFields(entry,'pairTransAngleZXZPsi','pairInvTransAngleZXZPsi');
entry=swapFields(entry,'pairTransAngleZXZTheta','pairInvTransAngleZXZTheta');

%Coordinates
entry=swapFields(entry,'pairCoordinateX1','pairCoordinateX2');
entry=swapFields(entry,'pairCoordinateY1','pairCoordinateY2');
entry=swapFields(entry,'pairCoordinateZ1','pairCoordinateZ2');

%Angles
entry=swapFields(entry,'pairAnglePhi1','pairAnglePhi2');
entry=swapFields(entry,'pairAnglePsi1','pairAnglePsi2');
entry=swapFields(entry,'pairAngleTheta1','pairAngleTheta2');

%AddInfo
entry=swapFields(entry,'pairClass1','pairClass2');
entry=swapFields(entry,'pairPsf1','pairPsf2');

function st=swapFields(st,field1,field2)

v1=getfield(st,field1);
v2=getfield(st,field2);

st=setfield(st,field1,v2);
st=setfield(st,field2,v1);

