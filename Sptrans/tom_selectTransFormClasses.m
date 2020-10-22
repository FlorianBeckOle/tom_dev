function [transListSel,selFolder]=tom_selectTransFormClasses(transList,selList,outputFolder)
%TOM_SELECTTRANSFORMCLASSES selects classes from transForm list
%
% transListSel= tom_selectTransFormClasses(transList,selList,outputFolder)
%
%PARAMETERS
%
%  INPUT
%  transList                transformation List
%  selList                    selecton List
%  outputFolder         folder for output
%
%  OUTPUT
%   transListSel        sleleted lists
%   selFolder            folder where selected lists have been written 
% 
%EXAMPLE
%
% %configure 
% selList(1).classNr=[1 2];
% selList(1).polyNr=-1;
%
% %run
%  tom_selectTransFormClasses(transList,selList,'run1')
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
    outputFolder='';
end;


    
if (ischar(transList))
    transList=tom_starread(transList);
end;
st=tom_extractData(transList);

if (ischar(selList))
    if (strcmp(selList,'Classes-Sep'))
        uClass=unique(st.label.pairClass(:));
        clear('selList');
        for i=1:length(uClass)
             selList(i).classNr=uClass(i);
             selList(i).polyNr=-1;
        end;
    end;
        
end;
 
 
 for i=1:length(selList)
     clNr=selList(i).classNr;
     polyNr=selList(i).polyNr;
     idx=filterList(st,clNr,polyNr);
     if (nargout>0)
         transListSel{i}=transList(idx);
     end;
     trFilt=transList(idx);
     trFilt(1).Header=transList(1).Header;
    selFolder{i}=genOutput(trFilt,selList(i),outputFolder);
 end;



function selFolder=genOutput(transList,selList,outputFolder)

clNr=selList.classNr;
polyNr=selList.polyNr;

if (polyNr(1)==-1)
    polyNrStr='';
else
    polyNrStr=strrep(num2str(polyNr),' ','+' );
    polyNrStr=['p' polyNrStr];
end;
selFolder=[outputFolder filesep 'c' strrep(num2str(clNr),' ','+') polyNrStr];
selFolder=strrep(selFolder,'++','+');

warning off;
mkdir(selFolder);
warning on;


tom_starwrite([selFolder filesep 'transList.star'],transList);


function idx=filterList(st,classNr,polyNr)
  
     if (isfield(st.label,'pairClass'))
         if (classNr==-1)
             idxC=1:length(st.label.pairClass);
         else
             allClasses=st.label.pairClass;
             idxC=find(ismember(allClasses,classNr));
         end;
     else
         idxC=1:size(st.p1.positions,1);
     end;
     
     if (isfield(st.label,'pairLable'))
         if (polyNr==-1)
             idxP=1:length(st.label.pairClass);
         else
             allPoly=st.label.pairLable;
             idxP=find(ismember(allPoly,polyNr));
         end;
     else
          idxP=1:size(st.p1.positions,1);
     end; 
     
 idx=intersect(idxP,idxC);

