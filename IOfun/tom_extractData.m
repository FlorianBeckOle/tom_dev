function st=tom_extractData(listFile,makePosUnique)
% tom_extractData reads positions,angles... from listFiles (.star)
%  
%     [angles,positions,shifts,cmbInd]=tom_extractData(listFile)
%  
%  PARAMETERS
%  
%    INPUT
%     listFile                  filename of the list
%     makePosUnique   (0) flag for removing duplicates
%
%    OUTPUT
%  
%      st                            structure containing list information 
%      
%
%     p1.angles              angles (nx3) in zxz
%     p1.positions          positions  (nx3)          
%     p1.shifts                shifts (nx3)
%     p1.classes             classes 
%     p1.tomoNames     tomoName (n) cell
%     p1.tomoIDs           tomoNames translated inot ids 
%     p1.psfs                  point spread for particles
%     p1.pixs                  pixelsize in Ang
%     p1.cmbInd            cmbInd index of  pairs only for pair files 
%
%
%  EXAMPLE
%       [angles,positions,shifts,classes]=tom_extractData('test.star');
%
%  REFERENCES
%  
%  NOTE:
%  
%
%  SEE ALSO
%      tom_xmippdocread,tom_spiderread
%  
%     created by FB 12/08/19
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom
% 

if (nargin<2)
    makePosUnique=0;
end;

if (iscell(listFile)==0)
    tmp=listFile;
    clear('listFile'); listFile{1}=tmp;
end;


for i=1:length(listFile)
    actList=listFile{i};
    st=extractData(actList,makePosUnique);
end;
st=updateTomoID(st);




function  st=updateTomoID(st)
    
tomoNames=st.p1.tomoName;
tomoIDs=zeros(length(tomoNames),1);
utomoName=unique(tomoNames);
for i=1:length(utomoName)
    idx= ismember(tomoNames,utomoName{i});
    tomoIDs(idx)=i;
end;
st.label.tomoID=tomoIDs;
    

function st=catData(st)


    angles=cat(1,angles,anglesAct);
    positions=cat(1,positions,positionsAct);
    shifts=cat(1,shifts,shiftsAct);
    classes=cat(1,classes,classesAct);
    tomoNames=cat(2,tomoNames,tomoNameAct);
    psfs=cat(2,psfs,psfsAct);
    pixs=cat(1,pixs,pixsAct);
    cmbInd=cat(1,cmbInd,cmbIndAct);



function st=extractData(listFile,makePosUnique)

 [list,type]=readList(listFile);

 if (strcmp(type,'pairStar'))
    st=extractFromPairStar(list);
 end;
 
 if (strcmp(type,'relionStar'))
     st=extractFromRelionStar(list);
 end;
 
% if (strcmp(ext,'.star'))
%     
%      
%      positions=zeros(length(star),3);
%      angles=zeros(length(star),3);
%      shifts=zeros(length(star),3);
%      classes=zeros(length(star),1);
%      pixs=zeros(length(star),1);
%      
%      if ( isempty(find(ismember(fieldnames(star),'pairTransVectX')))==0 )
%         
%      else
%          for i=1:length(star)
%              positions(i,:)=[star(i).rlnCoordinateX star(i).rlnCoordinateY star(i).rlnCoordinateZ];
%              [~,angles(i,:)] = tom_eulerconvert_xmipp(star(i).rlnAngleRot,star(i).rlnAngleTilt,star(i).rlnAnglePsi);
%              classes(i) = star(i).rlnClassNumber;
%              shifts(i,:) =[star(i).rlnOriginX star(i).rlnOriginY star(i).rlnOriginZ];
%              tomoName{i}=star(i).rlnMicrographName;
%              psfs{i}=star(i).rlnCtfImage;
%              pixs(i)=star(i).rlnDetectorPixelSize./star(i).rlnMagnification.*10000;
%          end;
%      end;
%     
% end;
%  
% 
% if (strcmp(ext,'.pair'))
%     clear('cmbInd');  
%     list=load(listFile);
%     cmbInd=zeros(size(list,1),3);
%     positions=zeros(size(list,1)*2,3);
%     angles=zeros(size(list,1)*2,3);
%     zz=0;
%      for i=1:size(list,1)
%         zz=zz+1;
%         positions(zz,:)=list(i,2:4);
%         angles(zz,:)=list(i,5:7);
%         zz=zz+1;
%         positions(zz,:)=list(i,9:11);
%         angles(zz,:)=list(i,12:14);
%         cmbInd(i,1)=zz-1;
%         cmbInd(i,2)=zz;
%         cmbInd(i,3)=list(i,8);
%     end;
% end;
% 
% if (makePosUnique)
%     [~,idx]=unique(positions,'rows');
%     positions=positions(idx,:);
%     angles=angles(idx,:);
%     if (isempty(classes)==0)
%         classes=classes(idx);
%     end;
%     if (isempty(shifts)==0)
%         shifts=shifts(idx,:);
%     end;
% end;


function [list,type]=readList(listFile)
 
if (ischar(listFile))
    [~,~,ext]=fileparts(listFile);
else
    ext='';
    if (isstruct(listFile))
        list=listFile;
        if ( isempty(find(ismember(fieldnames(list),'pairTransVectX'), 1))==0 )
            type='pairStar';
        end;
        if ( isempty(find(ismember(fieldnames(list),'rlnCoordinateX'), 1))==0 )
            type='relionStar';
        end;
    end;
end;

if (strcmp(ext,'.star'))
    if (ischar(listFile))
        list=tom_starread(listFile);
    else
        list=listFile;
    end
    if ( isempty(find(ismember(fieldnames(list),'pairTransVectX'), 1))==0 )
        type='pairStar';
    end;
     if ( isempty(find(ismember(fieldnames(list),'rlnCoordinateX'), 1))==0 )
        type='relionStar';
    end;
 end;

 if (strcmp(ext,'.pair'))
     list=load(listFile);
     type='pairTxt';
 end;



function st=extractFromRelionStar(star)

    
for i=1:length(star)
     st.p1.positions(i,:)=[star(i).rlnCoordinateX star(i).rlnCoordinateY star(i).rlnCoordinateZ];
     [~,st.p1.angles(i,:)] = tom_eulerconvert_xmipp(star(i).rlnAngleRot,star(i).rlnAngleTilt,star(i).rlnAnglePsi);
    st.p1.classes(i) = star(i).rlnClassNumber;
    st.p1.tomoName{i}=star(i).rlnMicrographName;
    st.p1.psfs{i}=star(i).rlnCtfImage;
    st.p1.pixs(i)=star(i).rlnDetectorPixelSize./star(i).rlnMagnification.*10000;
    st.label.tomoName{i}=star(i).rlnMicrographName;
end; 
    
    
  
 function st=extractFromPairStar(star)



cmbInd=zeros(length(star),7);
zz=0;
for i=1:length(star)
  
    st.p1.positions(i,:)=[star(i).pairCoordinateX1 star(i).pairCoordinateY1 star(i).pairCoordinateZ1];
    st.p1.angles(i,:) = [star(i).pairAnglePhi1 star(i).pairAnglePsi1 star(i).pairAngleTheta1];
    st.p1.classes(i) = star(i).pairClass1;
    st.p1.tomoName{i}= star(i).pairTomoName;
    st.p1.psfs{i}=star(i).pairPsf1;
    st.p1.pixs(i)=star(i).pairPixelSizeAng;
    st.p1.orgListIDX(i)=star(i).pairIDX1;
   
    st.p2.positions(i,:)=[star(i).pairCoordinateX2 star(i).pairCoordinateY2 star(i).pairCoordinateZ2];
    st.p2.angles(i,:) = [star(i).pairAnglePhi2 star(i).pairAnglePsi2 star(i).pairAngleTheta2];
    st.p2.classes(i) = star(i).pairClass2;
    st.p2.tomoName{i}= star(i).pairTomoName;
    st.p2.psfs{i}=star(i).pairPsf2;
    st.p2.pixs(i)=star(i).pairPixelSizeAng;
    st.p2.orgListIDX(i)=star(i).pairIDX2;
    
    st.label.pairClass(i)=star(i).pairClass;
    st.label.pairClassColour(i,:)=[str2num(strrep(star(i).pairClassColour,'-',' '))];
    st.label.pairLable(i,:)=star(i).pairLabel;
    st.label.tomoName{i}=star(i).pairTomoName;
    st.label.tomoID(i)=-1;
    st.label.p1p2TransVect(i,:)=[st.p2.positions(i,:)- st.p1.positions(i,:)];
    st.label.orgListName{i}=star(i).pairOriPartList;
    
    cmbInd(i,1)=star(i).pairIDX1;
    cmbInd(i,2)=star(i).pairIDX2;
    cmbInd(i,3)=star(i).pairClass;
    ColourVect=str2num(strrep(star(i).pairClassColour,'-',' '));
    cmbInd(i,4:6)=ColourVect;
    cmbInd(i,7)=star(i).pairLabel;
    zz=zz+1;
end;


st.path=cmbInd;

