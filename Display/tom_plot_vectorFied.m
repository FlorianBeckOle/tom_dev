function  tom_plot_vectorFied(posAng,tomoID,classNr,polyNr,onlySelected,scale,repVect,col,cmbInd,outputFolder)
%TOM_PLOT_VECTORFILED plots template matching/subtomoAvg results as vectorsField
%
%
%   tom_plot_vectorFied(posAng,repVect,scale,col,plotFlag)
%   
%
%PARAMETERS
%
%  INPUT
%   posAng               nx6 matrix of  positions=>posAng(:,1:3) and rotations=>posAng(:,4:6) 
%                               or relion .star file or pair.star file
%   tomoID                (-1) id of tomogram
%   classNr               (-1) or vetor of selected classes
%   polyNr                (-1) show only this polysomes with this poly ID
%   onlySelected      (1) show only selected classes use 0 to switch off
%   scale                   (50) scale of vectors  usually the size of the template in pixels
%   repVect              ([0 0 1]) matrix of vectors to represent the template
%   col                      ([0 0 1]) color for rep-vectors
%   cmbInd               ('') index which vectors are pairs format vect1 vect2 chainId
%   outputFolder      ('') folder for rendered images               
%
%  OUTPUT                   
%   -
%
%EXAMPLE
%   tom_plot_vectorFied('../../rawFromQiang/starAvg/118.star');
%   %render to disk
%   tom_plot_vectorFied('clusterAll/run0/allTransforms.star',1:18,-1,-1,1,50,[0 0 1;1 1 0],[0 0 1],'','out/');
%
%REFERENCES
%
%SEE ALSO
%   tom_sum_rotation
%
%   created by FB 04/05/07
%     
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
   tomoID=-1;
end;

if (nargin<3)
   classNr=-1;
end;

if (nargin<4)
   polyNr=-1;
end;

if (nargin<5)
    onlySelected=1;
end;

if (nargin<6)
    scale=50;
end;

if (nargin<7)
    repVect=[1 0 0];
end;

if (nargin<8)
    col=[0 0 1];
end;

if (nargin<9)
    cmbInd='';
end;

if (nargin<10)
    outputFolder='';
end;


if (isnumeric(posAng) || ischar(posAng) || isstruct(posAng) )
    plot_vectorFied(posAng,repVect,scale,col,cmbInd,classNr,polyNr,onlySelected,tomoID,outputFolder);
end;
    
if (iscell(posAng))   
    for i=1:length(posAng)
        posAngAct=posAng{i};
        figure;
        plot_vectorFied(posAngAct,repVect,scale,col,cmbInd,classNr,polyNr,onlySelected,tomoID,outputFolder);
    end;
end;



function plot_vectorFied(posAng,repVect,scale,col,cmbInd,classNr,polyNr,onlySelected,tomoID,outputFolder)


if (isnumeric(posAng))
    fTitle='';
    pos=posAng(:,1:3);
    angles=posAng(:,4:6);
    plotRepVects(pos,angles,repVect,scale,col,fTitle);
    return;
end;



if (ischar(posAng) || isstruct(posAng))
     
    fTitleList='';
    
     if (ischar(posAng))
        fTitleList=[posAng ' '];
    end;
       
    st=tom_extractData(posAng);
    allTomoID=st.label.tomoID;
    uTomoID=unique(allTomoID);
    allTomoLabel=st.label.tomoName;
    uTomoLabel=unique(allTomoLabel);
  
      if (tomoID>-1)
        uTomoID=tomoID;
      end;
    
    if (length(uTomoID)>1)
        disp('more than one tomogram in list');
        for i=1:length(uTomoID)
              disp([ 'tomoID: ' num2str(uTomoID(i))  ' tomoName: ' uTomoLabel{i} ] );
        end;
    end;
    
    if ((length(uTomoID)>5)  && (tomoID(1)==-1 ) && isempty(outputFolder) ) 
        disp(['warning found ' num2str(length(uTomoID)) ' tomograms reducing to 5' ]);
        disp(['use tomoID parameter to select specific tomograms']);
        uTomoID=uTomoID(1:5);
    end;
    
  plotClassZero=isempty(find(classNr==0))==0;
  
  
  wb=tom_progress(length(uTomoID),'rendering vector fields');
  
  if (isempty(outputFolder))
      for i=1:length(uTomoID)
          doRender(st,classNr,polyNr,uTomoID,outputFolder,i,onlySelected,fTitleList,uTomoLabel,repVect,scale,col,plotClassZero);
          wb.update();
      end;
  else
      parfor i=1:length(uTomoID)
          doRender(st,classNr,polyNr,uTomoID,outputFolder,i,onlySelected,fTitleList,uTomoLabel,repVect,scale,col,plotClassZero);
          wb.update();
      end;
  end;
  
  wb.close();
  
end;


function doRender(st,classNr,polyNr,uTomoID,outputFolder,i,onlySelected,fTitleList,uTomoLabel,repVect,scale,col,plotClassZero)
        
    
    if (isempty(outputFolder)==0)
        warning off; mkdir(outputFolder); warning on;
        h=figure;  h.Visible='off';
        set(h,'Name','vector representation');
    else
        if (i>1)
            figure;
        end;
    end;
    
    idx=filterList(st,classNr,polyNr,uTomoID(i));
    
    if (isempty(idx))
        close(gcf);
        return;
    end;
    
    if (onlySelected)
        idxRep=idx;
    else
        idxRep=1:size(st.p1.positions,1);
    end;
    
    pos=st.p1.positions(idxRep,:);
    angles=st.p1.angles(idxRep,:);
    
    fTitle=cat(2,fTitleList,uTomoLabel{uTomoID(i)});
    plotRepVects(pos,angles,repVect,scale,col,fTitle);
    
    if (isfield(st,'p2'))
        pos=st.p2.positions(idxRep,:);
        angles=st.p2.angles(idxRep,:);
        hold on;
        plotRepVects(pos,angles,repVect,scale,col,fTitle);
        hold on;
        plotPairs(st,idx,plotClassZero);
        hold off;
    end;
    fTitle='';
    if (isempty(outputFolder)==0)
        [~,fnameTmp]=fileparts(uTomoLabel{uTomoID(i)});
        h.Visible='on'; saveas(h,[outputFolder filesep fnameTmp],'fig'); close(h);
    end;
        


 function idx=filterList(st,classNr,polyNr,tomoID)
  
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
      
     if (isfield(st.label,'tomoID'))
         if (tomoID==-1)
             idxT=1:length(st.label.pairClass);
         else
             allTomo=st.label.tomoID;
             idxT=find(ismember(allTomo',tomoID));
         end;
     else
          idxT=1:size(st.p1.positions,1);
     end; 
     
     
 idx=intersect(idxP,idxC);
 idx=intersect(idx,idxT);
 
  
 
function plotPairs(st,idx,plotClassZero)

 allClasses=st.label.pairClass(idx);
 uClasses=unique(allClasses);
 allTrans=st.label.p1p2TransVect(idx,:);   
 allCol=st.label.pairClassColour(idx,:);
 allPos=st.p1.positions(idx,:);
 allLabel= st.label.pairLable(idx,:);
 
 for i=1:length(uClasses)
   
    if (uClasses(i)==0)
       if (plotClassZero==0)
            continue;
       end;
    end;
     
    tmpIdx=find(allClasses==uClasses(i));
   connVect=allTrans(tmpIdx,:);
   conCol=allCol(tmpIdx(1),:);
   conPos=allPos(tmpIdx,:);
   midPos= conPos+(connVect.*0.35);
   quiver3(conPos(:,1),conPos(:,2),conPos(:,3),connVect(:,1),connVect(:,2),connVect(:,3),0,'g--','color',conCol,'LineWidth',1,'MaxHeadSize',1,'AutoScale','off');
   axis image;
   
   for ii=1:length(tmpIdx)
        lablePoly{ii}=['c' num2str(uClasses(i))  ',p' num2str(allLabel(tmpIdx(ii)))  ]  ;
   end;
   
   text(midPos(:,1),midPos(:,2),midPos(:,3),lablePoly,'FontSize',8);   
   clear('lablePoly');
   
 end;
     
 


 function plotRepVects(pos,angles,repVect,scale,col,fTitle) 

 for i=1:size(repVect,1)
     vectTmp=repVect(i,:);
     vectsRot=tom_rotVectByAng(vectTmp,angles);
     quiver3(pos(:,1),pos(:,2),pos(:,3),vectsRot(:,1).*scale,vectsRot(:,2).*scale,vectsRot(:,3).*scale,0,'color',col,'AutoScale','off');
     axis image;
     if (isempty(fTitle)==0)
         title(fTitle);
     end;
     hold on;
 end;
 hold off;
     
