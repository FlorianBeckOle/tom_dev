function tom_genListFromTransForm(transList,outputFolder,posType,listFlav)
%TOM_GENLISTFROMTRANSFORM generate list form transform paris (e.g. relion .start file or .coords)
%
% transListSel= tom_selectTransFormClasses(transList,selList,outputFolder)
%
%PARAMETERS
%
%  INPUT
%  transList                      transformation List
%  outputFolder              folder for output
%  posType                     ('center') center of pair
%                                      'particle'  as center
%  listFlav                        ('rel') flavour for output 
%                                       rel for relion  
%
%
%  OUTPUT
%   tom_genListFromTransForm(transList,'class1','center');
%   tom_genListFromTransForm(transList,'class1','particle','rel');
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 10/24/19
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


if (nargin<3)
    posType='center';
end;

if (nargin<4)
    listFlav='rel';
end;



if (ischar(transList))
    transList=tom_starread(transList);
end;
st=tom_extractData(transList);

warning off; mkdir(outputFolder); warning on;


if (strcmp(posType,'center'))
    baseNameOut=[outputFolder filesep 'pairCenter']; 
    ext='.coords';
    writePairCenterCoords(st,baseNameOut,ext);
end;

if (strcmp(posType,'particle'))
   writeParticleCenterList(st,listFlav,outputFolder);
end;

% 
% fileNameOut=[outputFolder filesep 'pairCenter.star'];
% writePairCentersStar(st,fileNameOut);





function writePairCenterCoords(st,baseNameOut,ext)



allTomoNames=st.label.tomoName;
allTomoNamesU=unique(st.label.tomoName);


for tomoNr=1:length(allTomoNamesU)
    
   
    [~,name]=fileparts(allTomoNamesU{tomoNr});
    fileNamOut=[baseNameOut '_' name  ext];
    fileNamOutCmd=[baseNameOut '_' name  ext '.angles.zyz'];
    
    fid = fopen(fileNamOut,'w');
    fid2= fopen(fileNamOutCmd,'w');
    idx=find(ismember(allTomoNames,allTomoNamesU(tomoNr)));
    positionsP1= st.p1.positions(idx,:);
    positionsP2= st.p2.positions(idx,:);
    anglesP1=st.p1.angles(idx,:);
    anglesP2=st.p2.angles(idx,:);
    
    for i=1:size(positionsP1,1)
        p1=positionsP1(i,:); p2=positionsP2(i,:);
        pM=round([p1 + p2]./2);
        fprintf(fid,'%d %d %d\n',pM(1),pM(2),pM(3));
    
        a1=anglesP1(i,:);
        a2=anglesP2(i,:);
        aM=tom_average_rotations([a1; a2]);
        [~,aMZYZ]=tom_eulerconvert_xmipp(aM(1),aM(2),aM(3),'tom2xmipp'); 
        fprintf(fid2,'%d %d %d %f %f %f\n',pM(1),pM(2),pM(3),aMZYZ(1),aMZYZ(2),aMZYZ(3));
    end;
    fclose(fid);
    fclose(fid2);
end;



function writeParticleCenterList(st,listFlav,outFold)

starOrg=tom_starread(st.label.orgListName{1});

all_ind=cat(1,st.p1.orgListIDX,st.p2.orgListIDX);
all_ind=unique(all_ind);

starNew=starOrg(all_ind);
starNew(1).Header=starOrg(1).Header;
tom_starwrite([outFold filesep 'allPart.star'],starNew);




