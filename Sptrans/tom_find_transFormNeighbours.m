function pairList=tom_find_transFormNeighbours(pairList,outputName,nrStatOut)
%TOM_FIND_TRANSFORMNEIGHBOURS finds neighbouring transform classes
%
%   pairList=tom_find_transFormNeighbours(pairList,outputName)
%
%PARAMETERS
%
%  INPUT
%   pairList               pari star file                  
%   outputName       (opt.) name of the output pair star file
%   nrStatOut            (5) show the 5 most abundant neigh classes
%
%  OUTPUT
%   pairList          with neighbours
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
%   created by FB 01/24/19
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


if (nargin<2)
    outputName='';
end;

if (nargin<3)
    nrStatOut=10;
end;

if (ischar(pairList))
    pairList=tom_starread(pairList);
end;

allTomoId=[pairList(:).pairTomoID];
allTomoIdU=unique(allTomoId);


parfor i=1:length(allTomoIdU)
      idx{i}=find(allTomoId==allTomoIdU(i));
      listTomo{i}=pairList(idx{i});  
end;

for i=1:length(allTomoIdU)
  [listTomo{i},neighNpT{i},neighNMT{i}]=findNeighOneTomo(listTomo{i});
end;

neigh_N_puls=[];
neigh_N_minus=[];
for i=1:length(listTomo)
    pairList(idx{i})=listTomo{i};
    neigh_N_puls=cat(1,neigh_N_puls,neighNpT{i});
    neigh_N_minus=cat(1,neigh_N_minus,neighNMT{i});
end;

nCmb=cat(2,neigh_N_puls,neigh_N_minus);
nCmb=nCmb(find(ismember(nCmb,[0 0 0 0 0 0],'rows')==0),:);
[nCmbU,~,ic]=unique(nCmb,'rows');
[clCount]=hist(ic,unique(ic));
[clCount_Sort,clIdx_Sort]=sort(clCount,'descend');

if (nrStatOut>length(clCount_Sort))
    nrStatOut=length(clCount_Sort);
end;

disp('classn+1    ||    classn-1        abundance' )
for i=1:nrStatOut
    disp(['  ' num2str(nCmbU(clIdx_Sort(i),1:3)) '       ||       '  num2str(nCmbU(clIdx_Sort(i),4:6))  '           '  num2str(round(clCount_Sort(i)./length(pairList).*100,1)) '%' ]);
end;

if (isempty(outputName)==0)
    tom_starwrite(outputName,pairList);
end;


function tmp=vect2ClassStr(vect)

tmp=strrep(num2str(vect),'-1','');
tmp=strrep(tmp,' ',';');
tmp=strrep(tmp,';;',';');
if (isempty(tmp) || strcmp(tmp(1),';') )
    tmp=['-1;'];
end;
if (strcmp(tmp(end),';')==0)
   tmp=[tmp ';']; 
end;



function [classes,p_act]=find_neigh(all_idx,all_pos,all_class,act_Idx)

  p_act=[];
  tmpInd=find(all_idx==act_Idx);
    
    if (isempty(tmpInd))
        tmpInd=-1;
    end;
 
    
    if (tmpInd>0)
         p_act=all_pos(tmpInd,:);
         tmp=all_class(tmpInd);
         if (length(tmp)>3)
            tmp=tmp(1:3);
         end;
         if (length(tmp)<3)
                tt=ones(1,3).*-1;
                tt(1:length(tmp))=tmp;
                tmp=tt;
         end;
         classes=tmp;
    else
        classes=-1;
    end;

  classes(find(classes==0))=-1; 

  if (isempty(p_act)==0)
      p_act=p_act(1,:);
  end;
  
  

 function [pairList,neigh_N_puls,neigh_N_minus]=findNeighOneTomo(pairList)
  
 
debugFl=0;    
     
st=tom_extractData(pairList);
idxU=unique(cat(2,st.p1.orgListIDX,st.p2.orgListIDX));
neigh_N_puls=zeros(max(idxU),3);
neigh_N_minus=zeros(max(idxU),3);
for i=1:length(idxU)
  
     all_class=st.label.pairClass;
     act_idx=idxU(i);
     
     
     all_idx1=st.p1.orgListIDX;
     all_pos1=st.p1.positions;
     [classes1,pT1]=find_neigh(all_idx1,all_pos1,all_class,act_idx);
    
     all_idx2=st.p2.orgListIDX;
     all_pos2=st.p2.positions;
      [classes2,pT2]=find_neigh(all_idx2,all_pos2,all_class,act_idx);

      if (isempty(pT1)==0)
          pT=pT1;
      else
          pT=pT2;
      end;
      
    neigh_N_puls(idxU(i),:)=classes1;
    neigh_N_minus(idxU(i),:)=classes2;
   
    
    if (debugFl==1)
        hold on; plot3(pT(1),pT(2),pT(3),'ro'); hold off;
        hold on; text(pT(1),pT(2),pT(3),['p:' num2str(classes1) ' m: ' num2str(classes2)],'FontSize',8); hold off;
    end;
end;


idx1=[pairList(:).pairIDX1];
idx2=[pairList(:).pairIDX2];
for i=1:length(idxU)
     act_idx=idxU(i);
     
     tmpIdx1=find(idx1==act_idx);
     for ii=tmpIdx1
            pairList(ii).pairNeighPlus1=vect2ClassStr(neigh_N_puls(act_idx,:));
            pairList(ii).pairNeighMinus1=vect2ClassStr(neigh_N_minus(act_idx,:));
     end;
   
      tmpIdx2=find(idx2==act_idx);
      for ii=tmpIdx2
            pairList(ii).pairNeighPlus2=vect2ClassStr(neigh_N_puls(act_idx,:));
            pairList(ii).pairNeighMinus2=vect2ClassStr(neigh_N_minus(act_idx,:));
     end;
end;

