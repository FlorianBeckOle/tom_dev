function  starSt=tom_calcTransforms(posAng,tomoNames,maxDist,dmetric,outputName,verbose)
%TOM_CALCTRANSFORMS calcutes pair transformations for all combinations use
%pruneRad to reduce calc
%
%   transList=tom_calcTransforms(pos,pruneRad,dmetric,verbose)
%
%PARAMETERS
%
%  INPUT
%   posAng            positions and angles of particles/observation nx6 or List .star               
%   tomoNames     ('') cell with tomogramName to seperate tomograms
%   maxDist          (inf) max particle/observation distance 
%                          pairs above this distance will be discarted
%                          usually a good value is 2.5 particleRad 
%  dmetric           ('exact') at the moment only exact implemented
%  outputName    ('') name of pair output starFile 
%  verbose           (1) verbose flag use 0 for no output
%
%
%  OUTPUT
%   starSt          star struck  containing the transformations and input Data
%
%
%EXAMPLE
%  parpool('local',32);  
%  datMat=tom_calcTransforms([0 0 0;1 1 1;2 2 2],[0 0 0; 0 0 10; 10 20 30]);
%    
%
%REFERENCES
%
%SEE ALSO
%   tom_calcPairTransForm,
%
%   created by FB 16/08/19
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
    tomoID='';
end;

if (nargin<3)
    maxDist=inf;
end;

if (nargin<4)
    dmetric='exact';
end;
 
if (nargin<5)
    outputName='';
end;

if (nargin<6)
    verbose=1;
end;

oriPartList='noListGiven';
if (ischar(posAng))
    oriPartList=posAng;
    posAng=tom_starread(posAng);
end;

if (isstruct(posAng))
    if (isfield(posAng,'pairTransVectX'))
        starSt=posAng;
        return;
    end;
end;

if  (ischar(posAng))
   st=tom_extractData(posAng,0);
   % [angles,pos,~,classes,tomoNames,tomoID,psfs,pixs]=tom_extractData(posAng,0);
    oriPartList=posAng;
else
    st=tom_extractData(posAng,0);
%     if (isempty(tomoNames))
%         tomoNames=ones(size(pos,1),1);
%     end;
    
end;



uTomoId=unique(st.label.tomoID);
uTomoNames=unique(st.label.tomoName);
transList=[];
allTomoNames=[];
idxOffSet=0;
for i=1:length(uTomoId)
    idx=find(st.label.tomoID==uTomoId(i));
    posAct=st.p1.positions(idx,:);
    anglesAct=st.p1.angles(idx,:);
    transListAct=calcTransforms(posAct,anglesAct,maxDist,dmetric,uTomoId(i),verbose);
    transListAct(:,1:2)= transListAct(:,1:2)+idxOffSet;
    transList=cat(1,transList,transListAct);
    allTomoNames=cat(1,allTomoNames,repmat({uTomoNames{i}},size(posAct,1),1));
    idxOffSet=idxOffSet+size(posAct,1);
end;
starSt=genStarFile(transList,allTomoNames,st,maxDist,oriPartList,outputName);



function starSt=genStarFile(transList,allTomoNames,st,maxDist,oriPartList,outputName)

classes=st.p1.classes';
psfs=st.p1.psfs;
pixs=st.p1.pixs;


header. isLoop=1;
header.title='data_';
header.fieldNames={'_pairIDX1','_pairIDX2','_pairTomoID' ...
                                  '_pairTransVectX','_pairTransVectY','_pairTransVectZ'...
                                  '_pairTransAngleZXZPhi', '_pairTransAngleZXZPsi','_pairTransAngleZXZTheta'...
                                  '_pairInvTransVectX','_pairInvTransVectY','_pairInvTransVectZ'...
                                  '_pairInvTransAngleZXZPhi', '_pairInvTransAngleZXZPsi','_pairInvTransAngleZXZTheta'...
                                  '_pairLenTrans','_pairAngDist' ...
                                  '_pairCoordinateX1','_pairCoordinateY1','_pairCoordinateZ1'...
                                  '_pairAnglePhi1','_pairAnglePsi1','_pairAngleTheta1'...
                                  '_pairClass1','_pairPsf1' ...
                                   '_pairNeighPlus1','_pairNeighMinus1' ...
                                  '_pairCoordinateX2','_pairCoordinateY2','_pairCoordinateZ2'...
                                  '_pairAnglePhi2','_pairAnglePsi2','_pairAngleTheta2'...
                                  '_pairClass2','_pairPsf2' ...
                                   '_pairNeighPlus2','_pairNeighMinus2' ...
                                  '_pairTomoName','_pairPixelSizeAng' ...
                                  '_pairOriPartList' ...
                                  '_pairMaxDist','_pairClass','_pairClassColour','_pairLabel','_pairScore'};
 
                              
listTransTmp=mat2cell(transList(:,1:17),ones(size(transList,1),1),ones(17,1));
idxTmp=transList(:,1:2);
listPart1=mat2cell(transList(:,18:23),ones(size(transList,1),1),ones(6,1));
listPart2=mat2cell(transList(:,24:29),ones(size(transList,1),1),ones(6,1));
classesPart1=mat2cell(classes(idxTmp(:,1)),ones(size(transList,1),1),ones(1,1));
classesPart2=mat2cell(classes(idxTmp(:,2)),ones(size(transList,1),1),ones(1,1));
psfsPart1=psfs(idxTmp(:,1));
psfsPart2=psfs(idxTmp(:,2));
neighPMPart{1}='-1;-1';
neighPMPart{2}='-1;-1';
neighPMPart=repmat(neighPMPart,size(transList,1),1);



tomoName12=allTomoNames(idxTmp(:,1));

addInfo{1}=pixs(1);
addInfo{2}=oriPartList;    
addInfo{3}=maxDist;
addInfo{4}=-1;
addInfo{5}='0-0-0';          
addInfo{6}=-1;
addInfo{7}=-1;          
  

addInfo=repmat(addInfo,size(transList,1),1);           
 
dataMat=listTransTmp;
dataMat=cat(2,dataMat,listPart1,classesPart1,psfsPart1',neighPMPart);
dataMat=cat(2,dataMat,listPart2,classesPart2,psfsPart2',neighPMPart);
dataMat=cat(2,dataMat,tomoName12);
dataMat=cat(2,dataMat,addInfo);

if (isempty(outputName)==0)
    tom_starwrite(outputName,dataMat,header);
end;

%TOD replace reading by in memory transform!!
%isStrVect=zeros(41,1);
%starSt=Raw2Struct(dataMat,isStrVect,header.fieldNames);

starSt=tom_starread(outputName);




function st=Raw2Struct(dataRaw,isStrVect,fieldnames)

structCommndString='st=struct(';
for i=1:length(fieldnames)
    if (isStrVect(i))
        structCommndString=[structCommndString '''' fieldnames{i}(2:end) '''' ',dataRaw{' num2str(i) '},' ];
    else
        structCommndString=[structCommndString '''' fieldnames{i}(2:end) '''' ', num2cell(dataRaw{' num2str(i) '}),' ]; 
    end;
end;
structCommndString=structCommndString(1:end-1);
structCommndString=[structCommndString ');'];
eval(structCommndString);


function  transList=calcTransforms(pos,angles,pruneRad,dmetric,tomoID,verbose)

jobList=zeros(size(pos,1).*size(pos,1),2);
zz=1;

for i=1:size(pos,1)
    for ii=(i+1):size(pos,1)
        pos1=pos(i,:);
        pos2=pos(ii,:);
        if (norm(pos1-pos2)>pruneRad)
            continue;
        end;
        jobList(zz,1)=i;
        jobList(zz,2)=ii;
        zz=zz+1;
    end;
end;
jobListN=jobList(1:zz-1,:);
jobList=jobListN;



if (verbose==1)
    waitbar=tom_progress(size(jobList,1)/10,['Tomo: ' num2str(tomoID) 'Calculating transforms of ' num2str(size(jobList,1)) ' pairs' ]);
end;

transList=zeros(size(jobList,1),29,'single');
%for i=1:size(jobList,1)
parfor i=1:size(jobList,1)
    icmb=jobList(i,:);
    pos1=pos(icmb(1),:);
    pos2=pos(icmb(2),:);
    ang1=angles(icmb(1),:);
    ang2=angles(icmb(2),:);
    
    [tr,ang,dtr,dang]=tom_calcPairTransForm(pos1,ang1,pos2,ang2,dmetric);
    [Inv_tr,Inv_ang]=tom_calcPairTransForm(pos2,ang2,pos1,ang1,dmetric);
    transList(i,:)=[icmb tomoID tr ang Inv_tr Inv_ang dtr dang pos1 ang1 pos2 ang2];
   
    if ((mod(i,10)==1 || i==size(jobList,1)) && verbose==1 )
       waitbar.update();
    end;
end;

if (verbose==1)
    waitbar.close();
    clear('watibar');
    disp(' ');
end;


