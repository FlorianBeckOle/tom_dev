function pairListAlg=tom_find_connectedTransforms(pairList,outputName,branchDepth)
%TOM_FIND_CONNECTEDTANSFORMS finds connedted transform 
%
%   pairs=tom_find_connectedTransforms(pairList,outputName)
%
%PARAMETERS
%
%  INPUT
%   pairList               pari star file                  
%   outputName       (opt.) name of the output pair star file
%   branchDepth     (1) (opt.) Depth for branch tracing
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
%   created by FB 01/24/06
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
    branchDepth=1;
end;

%parse inputs
if (ischar(pairList))
    pairList=tom_starread(pairList);
end;


allClasses=[pairList(:).pairClass];
allClassesU=unique(allClasses);

allTomos=[pairList(:).pairTomoID];
allTomosU=unique(allTomos);

br=zeros(length(allClassesU),1);
pairListAlg=pairList;
waitbar=tom_progress(length(allClassesU),'tracing polysomes');
for i=1:length(allClassesU)
     if (allClassesU(i)==0)
        continue;
    end;
    idx1=find(allClasses==allClassesU(i));
    
    offset_PolyID=0;
    for ii=1:length(allTomosU)
        idx2=find(allTomos==allTomosU(ii));
        idx=intersect(idx1,idx2);
        %disp(['processing '  num2str(allClassesU(i))  ' ' num2str(allTomosU(ii))  ' len ' num2str(length(idx)) ]);
        if ( ~isempty(idx) )
            [pairListAlg(idx),br(allClassesU(i)),offset_PolyID]=findConn(pairList(idx),branchDepth,offset_PolyID);
        end;
    end;
    waitbar.update();
end;
close(waitbar);

if (isempty(find(br))==0 )
    disp(['warning found branches in class nr: ' num2str(find(br)') ]);
    disp(['==>make smaller classes']);
end;

if (nargin>1)
    tom_starwrite(outputName,pairListAlg);
end;


function [pairList,branchFound,offset_PolyID]=findConn(pairList,brDepth,offset_PolyID)

ind1=[pairList(:).pairIDX1];
ind2=[pairList(:).pairIDX2];

if  ( length(pairList)>1 )
    
    idxboth=unique(cat(1,ind1,ind2));
    cmbInd=[ind1(:) ind2(:) zeros(size(ind2))'];
    
    branchNr=length(cmbInd(:,1))- length(unique(cmbInd(:,1)))+1;
    if (branchNr>1)
        branchFound=1;
        branchDepth=brDepth;
    else
        branchFound=0;
        branchDepth=1;
    end;
    
    allPathN={};
    for br=1:branchDepth
        for i=1:size(cmbInd,1)
            zz=i;
            try
                [tmpPath]=searchPathForward(cmbInd,zz,br);
            catch
                disp('  ');
                disp('error search Path');
                disp('  ');
                continue;
            end;
            
            allPathN=uniquePathAdd(allPathN,tmpPath);
            clear('tmpPath');
        end;
    end;
    
     for i=1:length(allPathN)
        tmp=allPathN{i};
        cmbInd(tmp(:,3),3)=i+offset_PolyID;
    end;
    
else
    cmbInd(1,1:3)=[1 1 (offset_PolyID+1)];
    branchFound=0;
end;


uLabels=unique(cmbInd(:,3));
for i=1:length(uLabels)
    idxTmp=find(uLabels(i)==cmbInd(:,3));
    for ii=1:length(idxTmp)
         pairList(idxTmp(ii)).pairLabel=uLabels(i);
    end;
end;
offset_PolyID=max(uLabels);

function [tmpPath,hasBranch]=searchPathForward(cmbInd,zz,branchNr)


zzPath=1;
tmpPath=[cmbInd(zz,1:2) zz];
tmpPath=uint64(tmpPath);
circRepeat=0;

hasBranch=0;
for ii=1:1000
    idx2Search=cmbInd(zz,1);
    idxT1=find(cmbInd(:,2)==idx2Search);
    
    if (length(idxT1)>1)
        idxT1=idxT1(branchNr);
        hasBranch=0;
    end;
    
    zz=idxT1;
    if (ismember([cmbInd(zz,1:2) zz],tmpPath,'rows'))
           circRepeat=1;
    end;
    
    if (isempty(idxT1) || circRepeat==1)
        break;
    end;
    zzPath=zzPath+1;
    tmpPath(zzPath,:)=[cmbInd(zz,1:2) zz];
    
end;

if (ii==1000)
    disp('max trials reached');
end;


function allPath=uniquePathAdd(allPath,newPath)

if (isempty(newPath))
    return;
end;

if (isempty(allPath))
    [~,idx]=sort(newPath(:,3)');
    allPath{1}=newPath(idx,:);
    return;
end;

memAN=zeros(length(allPath),1,'uint8');
memNA=zeros(length(allPath),1,'uint8');
memBrach=zeros(length(allPath),1,'uint8');

newPathP=sort(newPath(:,3));
for i=1:length(allPath)
    actPath=allPath{i};
    actPath=actPath(:,3);
   
    interSNPathactPath=fastIntersect(newPathP,actPath);
    memAN(i)=length(fastIntersect(actPath,newPathP))==length(newPathP);
    memNA(i)=length(interSNPathactPath)==length(actPath);
    memBrach(i)=0;
    
    if ( isempty(interSNPathactPath)==0 &&  length(interSNPathactPath)<max([length(actPath) length(newPath)]) )
        memBrach(i)=1;
    end;
end;

if  (isempty(find(memBrach))==0 ) 
       ind=find(memBrach);
       for iii=1:length(ind)
            newPathUnion=union(newPath,allPath{ind(iii)},'rows');
       end;
        [~,idx]=sort(newPathUnion(:,3)'); 
       
        for iii=1:length(ind)
            allPath{ind(iii)}=newPathUnion(idx,:);
        end;
       return
end;   
  
 if ((isempty(find(memAN))) && (isempty(find(memNA))) )
        [~,idx]=sort(newPath(:,3)');
        allPath{length(allPath)+1}=newPath(idx,:);
 end;
 
 if  (isempty(find(memNA))==0 ) 
       [~,idx]=sort(newPath(:,3)'); 
       allPath{find(memNA)}=newPath(idx,:);
 end; 


    function C= fastIntersect(A,B)
        
        if ~isempty(A)&&~isempty(B)
            P = zeros(1, max(max(A),max(B)) ) ;
            P(A) = 1;
            C = B(logical(P(B)));
        else
            C = [];
        end
        
 