function stat=tom_analyseTransFromPopulation(pairList,classNr,repVolume,outputFolder,verbose)
%TOM_ANALYSETRANSFORMPOPULATION anlyses a set of given transformations
%
%   stat=tom_analyseTransFromPopulation(posAngList,classNr,template,outputFolder,verbose)
%
%PARAMETERS
%
%  INPUT
%  pairList             list of transform pairs
%  classNr             (-1) class number to process use -1 for all classes 
%  repVolume       ('')  volume to represent transformation  
%  outputFolder    ('') path to output folderf
%  verbose            1 verbose flag
%
%  OUTPUT
%   stat           struct containig information about 
%
%EXAMPLE
%
%  stat=tom_analyseTransFromPopulation(pairList,-1,'','out');
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
%

if (nargin<2)
    classNr=-1;
end;

if (nargin<3)
    repVolume='';
end;

if (nargin<4)
    outputFolder='';
end;

if (nargin<5)
    verbose=1;
end;

if (isstruct(pairList)==0)
    pairList=tom_starread(pairList);
end;

allClasses=[pairList(:).pairClass];
allClassesU=unique(allClasses);

%pairListAlg=pairList;
for i=1:length(allClassesU)
    idx=find(allClasses==allClassesU(i));
   stat(i)=analysePopulation(pairList(idx));
end;
 stat=sortStat(stat);

writeOutputStar(stat,outputFolder);

genOutput(stat,2); 
% if (nargin>1)
%     tom_starwrite(outputName,pairListAlg);
% end;

function stat=sortStat(stat)

polyFact=([stat(:).numPolybg3]+[stat(:).numPolybg5])*10000+[stat(:).num];
[~,idx]=sort(polyFact,'descend');
stat=stat(idx);


function writeOutputStar(stat,outputFolder)

if (isempty(outputFolder)==0)
    Header.title='data_';
    Header.isLoop=1;
    Header.fieldNames=fieldnames(stat);
    stat(1).Header=Header;
    tom_starwrite([outputFolder filesep 'statPoly.star'],stat);
end;


function genOutput(stat,minClassMembers)

if (length(stat)>20)
    idx=find([stat.num]>2);
    stat=stat(idx);
end;

zz=0;


zz=zz+1;
h.txt{zz}='classNr  ';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['%' num2str(length(h.txt{zz})) 'd ' ];
data{:,zz}=[stat(:).classNr];

zz=zz+1;
h.txt{zz}='num   ';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['%' num2str(length(h.txt{zz})) 'd ' ];
data{:,zz}=[stat(:).num];

zz=zz+1;
h.txt{zz}='stdTrans';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['% ' num2str(length(h.txt{zz})+6) '.1f  '];
data{:,zz}=[stat(:).stdTransVect];

zz=zz+1;
h.txt{zz}='  stdAng';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['% ' num2str(length(h.txt{zz})+1) '.1f  '];
data{:,zz}=[stat(:).stdTransAng];


zz=zz+1;
h.txt{zz}=' nrPolybg5';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['%' num2str(length(h.txt{zz})+3) 'd ' ];
data{:,zz}=[stat(:).numPolybg5];

zz=zz+1;
h.txt{zz}=' nrPolybg3';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['%' num2str(length(h.txt{zz})+5) 'd ' ];
data{:,zz}=[stat(:).numPolybg3];

zz=zz+1;
h.txt{zz}=' nrPolyMax';
h.format{zz}=['%' num2str(length(h.txt{zz})) 's '];
d.format{zz}=['%' num2str(length(h.txt{zz})+3) 'd ' ];
data{:,zz}=[stat(:).numPolyMax];

%h.txt{3}='stdAng';
%h.txt{4}='stdAng';



fHead='';
fData='';
for i=1:length(h.txt)
    fHead=cat(2,fHead,h.format{i});
    fData=cat(2,fData,d.format{i});
end;
fHead=cat(2,fHead,'\n');
fData=cat(2,fData,'\n');


fprintf(fHead,h.txt{1},h.txt{2},h.txt{3},h.txt{4},h.txt{5},h.txt{6},h.txt{7});
for i=1:size(data{1},2)
    fprintf(fData,data{1}(i),data{2}(i),data{3}(i),data{4}(i),data{5}(i),data{6}(i),data{7}(i));
end;

if (length(stat)>20)
    disp('only classes with more than 2 transformations showed');
end;


disp('')

function stat=analysePopulation(pairList)


classNr=pairList(1).pairClass;

stat.classNr=classNr;
stat.num=length(pairList);

vects=[pairList(:).pairTransVectX; pairList(:).pairTransVectY; pairList(:).pairTransVectZ]';
if (size(vects,1)>1)
    meanV=mean(vects);
else
    meanV=vects;
end;

diffV=vects-repmat(meanV,size(vects,1),1);
 for i=1:size(diffV,1)
      lendiffV(i)=norm(diffV(i,:));
end;
stdTransVect=std(lendiffV);
stat.meanTransVectX=meanV(1);
stat.meanTransVectY=meanV(2);
stat.meanTransVectZ=meanV(3);
stat.stdTransVect=stdTransVect;

angs=[pairList(:).pairTransAngleZXZPhi; pairList(:).pairTransAngleZXZPsi; pairList(:).pairTransAngleZXZTheta]';
meanAng=tom_average_rotations(angs);
for i=1:size(angs,1)
        lendiffAng(i)=tom_angular_distance(angs(i,:),meanAng);
end;
stdTransAng=std(lendiffAng);
%stat.meanTransAng=meanAng;
stat.meanTransAngPhi=meanAng(1);
stat.meanTransAngPsi=meanAng(2);
stat.meanTransAngTheta=meanAng(3);
stat.stdTransAng=stdTransAng;


allLabel=[pairList(:).pairLabel];
allLabelU=unique(allLabel);

if (allLabelU(1)==-1)
   stat.numPolybg5=-1;
   stat.numPolybg3=-1;
   stat.numPolyMax=-1; 
else
    for i=1:length(allLabelU)
        allNum(i)=length(find(allLabel==allLabelU(i)));
    end;
   stat.numPolybg5=length(find(allNum>5));
   stat.numPolybg3=length(find(allNum>3));
   stat.numPolyMax=max(allNum); 
end




