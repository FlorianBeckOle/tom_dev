function call_clusterTransForms()



conf=tom_clusterTransForms();

%% computing
numOfProc=32; %use -1 to switch off


%%
doPass=[1]; 

%% initial clustering --pass 0

%io
conf.io.posAngList='cmbAll.star';
conf.io.projcetFolder='clusterAll';
conf.io.classificationRun='run0';

%transFormations
conf.transForm.pixS=3.42; %in Ang
conf.transForm.maxDist=342; %in Ang

%classify
conf.classify.clustThr='auto';
%conf.classify.clustThr=0.1;

%select classes for output
conf.sel='Classes-Sep'; %write all list for all classes

%visualisation
conf.vis.vectField.showTomo=-1;
conf.vis.vectField.showClass=-1;
conf.vis.vectField.onlySelected=1;
conf.vis.vectField.polyNr=-1;
conf.vis.vectField.repVect=[0 1 0];

zz=1;
allConf(zz)=conf;

%% fine clustering --pass 1 

%io
conf.io.posAngList='clusterAll/run0/classes/c1/transList.star';
conf.io.projcetFolder='clusterAll';
conf.io.classificationRun='run1';

%transFormations
conf.transForm.pixS=3.42; %in Ang
conf.transForm.maxDist=342; %in Ang

%classify
conf.classify.clustThr=0.18; %needs to be opt by the user 

%select classes for output
conf.sel='Classes-Sep'; %write all list for all classes

%visualisation
conf.vis.vectField.showTomo=-1;
conf.vis.vectField.showClass=-1;
conf.vis.vectField.onlySelected=1;
conf.vis.vectField.polyNr=-1;
conf.vis.vectField.repVect=[0 1 0];

zz=zz+1;
allConf(zz)=conf;

%% fine clustering --pass 2 

%io
conf.io.posAngList='clusterAll/run1/classes/c1/transList.star';
conf.io.projcetFolder='clusterAll';
conf.io.classificationRun='run2';

%transFormations
conf.transForm.pixS=3.42; %in Ang
conf.transForm.maxDist=342; %in Ang

%classify
conf.classify.clustThr=0.26; %needs to be opt by the user 

%select classes for output
conf.sel='Classes-Sep'; %write all list for all classes

%visualisation
conf.vis.vectField.showTomo=-1;
conf.vis.vectField.showClass=-1;
conf.vis.vectField.onlySelected=1;
conf.vis.vectField.polyNr=-1;
conf.vis.vectField.repVect=[0 1 0];

zz=zz+1;
allConf(zz)=conf;


%%


%% Code

p = gcp('nocreate');
if (isempty(p))
    disp('no pool created. starting parpool');
    parpool('local',numOfProc);
else
    disp('found parpool');
end;

for i=1:length(allConf)
    if (isempty(find(ismember(doPass,i-1)))==0) 
        tom_clusterTransForms(allConf(i));
    end;
end;

%%