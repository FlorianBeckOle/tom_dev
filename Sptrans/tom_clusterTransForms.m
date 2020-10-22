function pSt=tom_clusterTransForms(pSt)
%TOM_CLUSTERTRANSFORMS clusters transformations for given ang/position list
%
%   tom_clusterTransForms(paramst)
%
%PARAMETERS
%
%  INPUT
%  pSt             parameter structure
%  
%  OUTPUT
%   pSt           default parameter struct
%
%EXAMPLE
%
% %call use script recommended
% tom_clusterTransForms('script');
%
% %configure manual
% pSt=tom_clusterTransForms();
% pSt.io.posAngList='myList.star';
% pSt.classify.clustThr=0.0016;
%
% %run
%  tom_clusterTransForms(pSt);
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
%


%% main

if (nargin<1)
    disp('generating default paramters')
    pSt=genDefautStruct();
    return;
end;

if (ischar(pSt))
    pScript=which('call_clusterTransForms');
    disp('copy script to local folder');
    copyfile(pScript,'.'); 
    pSt='';
    return;
end;

pSt=initStructs(pSt);

transList=calcTransForms(pSt);
[transList,pSt]=groupTransForms(pSt,transList);

tom_analyseTransFromPopulation(transList,-1,'',[pSt.io.classifyFold filesep 'stat' filesep]);

tom_find_transFormNeighbours(transList,[pSt.io.classifyFold filesep 'allTransforms.star']);
[transListSel,selFolds]=tom_selectTransFormClasses(transList,pSt.sel, [pSt.io.classifyFold filesep 'classes']);

genOutputLists(pSt,transListSel,selFolds);

visResult(pSt,transList);


%% functions
function pSt=initStructs(pSt)
%[~,listBase]=fileparts(pSt.io.posAngList);
%pSt.io.projcetFolder=['cluster-' listBase];
pSt.io.classifyFold=[pSt.io.projcetFolder filesep pSt.io.classificationRun];
genOutpuFolder(pSt);



function pSt=genDefautStruct(pSt)

%% io
pSt.io.posAngList='list.star';
[~,listBase]=fileparts(pSt.io.posAngList);
%if (isempty(pSt.io.projcetFolde))
    pSt.io.projcetFolder=['cluster-' listBase];
%end;

pSt.io.classificationRun='run0';

%% trans
pSt.transForm.pixS=3.42;
pSt.transForm.maxDist=342;

%% classify
pSt.classify.clustThr=0.0016;

%% Select
pSt.sel(1).classNr=-1;
pSt.sel(1).polyNr=-1;


%% vis
pSt.vis.vectField.showTomo=-1;
pSt.vis.vectField.showClass=1:10000;
pSt.vis.vectField.onlySelected=1;
pSt.vis.vectField.polyNr=-1;
pSt.vis.vectField.repVect=[0 1 0];


function genOutputLists(pSt,transListSel,outputFoldSel)

wb=tom_progress(length(transListSel),'generating selection lists');
parfor i=1:length(transListSel)
    transListTmp=transListSel{i};
    outputFoldCenter=[outputFoldSel{i} filesep 'pairCenter/'] ;
    tom_genListFromTransForm(transListTmp,outputFoldCenter,'center');
    outputFoldCenter=[outputFoldSel{i} filesep 'particleCenter/'] ;
    tom_genListFromTransForm(transListTmp,outputFoldCenter,'particle');
    wb.update();
end;
wb.close();


function transList=calcTransForms(pSt)

maxDistInPix=pSt.transForm.maxDist./pSt.transForm.pixS;
transFormFile=[pSt.io.classifyFold filesep  'allTransforms.star'];

if (exist(transFormFile,'file') )
    disp(['load distances from: ' transFormFile]);
    transList=tom_starread(transFormFile);
else   
    transList=tom_calcTransforms(pSt.io.posAngList,'',maxDistInPix,'exact',transFormFile);
end;



function [transList,pSt]=groupTransForms(pSt,transList)

outputFold=[pSt.io.classifyFold filesep 'scores'];

treeFile=[outputFold filesep 'tree.mrc' ];

if (exist(treeFile,'file') )
    ll=tom_mrcread(treeFile);
    ll=ll.Value;
else
    ll=calcLinkage(transList,outputFold);
    tom_mrcwrite(ll,'name',treeFile);
end;

[clusters,~,~,thres]=tom_dendrogram(ll,pSt.classify.clustThr,size(transList,1),0); 
pSt.classify.clustThrAct=thres;


if (isempty(clusters))
    disp('warnning no classes check threshold');
    dendrogram(ll);
else
    for i=1:length(clusters)
        idx=clusters(i).members;
        class=clusters(i).id;
        if (idx==-1)
            continue;
        end;
        colour=[num2str(clusters(i).color(1)) '-' num2str(clusters(i).color(2)) '-' num2str(clusters(i).color(3))];
        for ii=1:length(idx)
            transList(idx(ii)).pairClass=class;
            transList(idx(ii)).pairClassColour=colour;
        end;
    end;
    transList=tom_align_transformDirection(transList);
    transList=tom_find_connectedTransforms(transList,[pSt.io.classifyFold filesep 'allTransforms.star']);
    
end;



function genOutpuFolder(pSt)

projcetFolder=pSt.io.projcetFolder;
classificationFolder=[projcetFolder filesep  pSt.io.classificationRun];

if (exist(projcetFolder,'dir')==0)
    mkdir(projcetFolder);
end

if (exist(classificationFolder,'dir')==0)
    mkdir(classificationFolder);
    mkdir([classificationFolder filesep 'scores']);
end

warning off;  mkdir([pSt.io.classifyFold filesep 'vis']); warning on;
warning off;  mkdir([pSt.io.classifyFold filesep 'vis/vectfileds/' ]); warning on;
warning off;  mkdir([pSt.io.classifyFold filesep 'vis/clustering/' ]); warning on;
warning off;  mkdir([pSt.io.classifyFold filesep 'stat']); warning on;

function visResult(pSt,transList)


disp(' ');
disp('rendering figures');

vectVisP=pSt.vis.vectField;
tom_plot_vectorFied(transList,vectVisP.showTomo,vectVisP.showClass,vectVisP.polyNr,vectVisP.onlySelected,50,vectVisP.repVect,[0 0 1],'',[pSt.io.classifyFold filesep 'vis/vectfileds/']);

nrTrans=length(transList);
treeFile=[pSt.io.classifyFold filesep 'scores/tree.mrc'];
thres=pSt.classify.clustThrAct;

parfor i=1:2
    if  (i==1)
          dspTree(pSt,nrTrans,treeFile);
      end;
    
    if  (i==2)
          dspLinkage(pSt,treeFile,thres);
    end;
end;

disp('rendering figures done');


function dspTree(pSt,nrTrans,treeFile)
disp(' ');
h=figure; h.Visible='off'; set(h,'Name','clustering');
tom_dendrogram(treeFile,pSt.classify.clustThr,nrTrans);
h.Visible='on'; saveas(h,[pSt.io.classifyFold filesep 'vis/clustering/tree'],'fig'); close(h);

function dspLinkage(pSt,treeFile,thres)

h=figure; h.Visible='off';  set(h,'Name','link-levels');
tree=tom_mrcread(treeFile);
tree=tree.Value;
plot(sort(tree(:,3),'descend' ));
hold on; plot(sort(tree(:,3),'descend' ),'ro' );
plot(ones(length(tree(:,3)),1).*thres,'k--');
legend({'link-level', 'link-level','threshold'});
h.Visible='on'; saveas(h,[pSt.io.classifyFold filesep 'vis/clustering/linkLevel'],'fig'); close(h);



function [ll]=calcLinkage(transList,preCalcFold)


transVect=[transList(:).pairTransVectX; transList(:).pairTransVectY; transList(:).pairTransVectZ]';
transVectInv=[transList(:).pairInvTransVectX; transList(:).pairInvTransVectY; transList(:).pairInvTransVectZ]';
distsVect=tom_pdist(transVect,'euc',transVectInv);
tom_mrcwrite(distsVect,'name',[preCalcFold filesep 'distsVectors.mrc']);


transAngVect=[transList(:).pairTransAngleZXZPhi; transList(:).pairTransAngleZXZPsi; transList(:).pairTransAngleZXZTheta]';
transAngVectInv=[transList(:).pairInvTransAngleZXZPhi; transList(:).pairInvTransAngleZXZPsi; transList(:).pairInvTransAngleZXZTheta]';
distsAng=tom_pdist(transAngVect,'ang',transAngVectInv);
tom_mrcwrite(distsAng,'name',[preCalcFold filesep 'distsAngles.mrc']);


if (std(distsVect)>0)
    distsVectNorm=tom_norm(tom_norm(distsVect,'mean0+1std'),1)+1;
else
    distsVectNorm=tom_norm(distsVect-mean(distsVect(:)),1)+1;
end;

if (std(distsAng)>0)
    distsAngNorm=tom_norm(tom_norm(distsAng,'mean0+1std'),1)+1;
else
     distsAngNorm=tom_norm(distsAng-mean(distsAng(:)),1)+1;
end;


distsCN=distsVectNorm.*distsAngNorm;
if (std(distsCN)>0)
    distsCN=tom_norm(tom_norm(distsCN,'mean0+1std'),1);
else
    distsCN=tom_norm(distsCN-mean(distsCN(:)),1)+1;
end;

tom_mrcwrite(distsCN,'name',[preCalcFold filesep 'distsCombined.mrc']);


ll=linkage(distsCN,'average');
tom_mrcwrite(ll,'name',[preCalcFold filesep 'tree.mrc']);



