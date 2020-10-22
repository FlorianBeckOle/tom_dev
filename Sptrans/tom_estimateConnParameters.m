function tom_estimateConnParameters()

%% io

io.trainList{1}='proj102/Classification/run1/allTransforms.star';
io.useClass(1)=10;
io.outPutTrainFile='output/train/train3Poly.mat';
io.outPutTransFormListAlg='proj102/Classification/run1/allTransformsAlg.star';


%% Clean Parameters
transPairsParam.clean.do=1;
transPairsParam.clean.nStd=2.5;
transPairsParam.clean.iter=15;

angPairsParam.clean.do=1;
angPairsParam.clean.nStd=2.5;
angPairsParam.clean.iter=15;

%% Distribution Parameters
transDist.type='normal';
transDist.prenorm='1std';
transDist.nStd=3.5;

angDist.type='normal';
angDist.prenorm='1std';
angDist.nStd=3.5;


%% Vis
vis.plotCleanUp=1;
vis.plotDistributions=1;
vis.plotTrainInput=1;

%% Code
pairs=readTrainData(io,vis);
pairs=alignToOneDirection(pairs);
pairsClean=cleanTrainData(pairs,transPairsParam,angPairsParam,vis);

disp('  estimating Distributions');
trainSt.pdfVect=genVectPdf(pairsClean,transDist,vis);
trainSt.pdfAng=genAngPdf(pairsClean,angDist,vis);

disp('  saving Distributions');
disp(['saving to: ' io.outPutTrainFile]);
save(io.outPutTrainFile,'trainSt');

analyseDists(trainSt);

function analyseDists(trainSt)

disp(' ');
disp(' score stat:');

if (strcmp(trainSt.pdfAng.prenorm ,'1std'))
    tmpStd=trainSt.pdfAng.std;
else
    tmpStd=1;
end;
scoresAng=pdf(trainSt.pdfAng.pdf,trainSt.pdfAng.data./tmpStd);

if (strcmp(trainSt.pdfVect.prenorm ,'1std'))
    tmpStd=trainSt.pdfVect.std;
else
    tmpStd=1;
end;
scoresVect=pdf(trainSt.pdfVect.pdf,trainSt.pdfVect.data./trainSt.pdfVect.std).*10;

trainScores=scoresAng.*scoresVect;
tom_dev(trainScores);
disp(['good thr. to start: ' num2str(min(trainScores)) ]);



function pdfST=genAngPdf(pairs,angDist,vis)

if ( isempty(find(ismember(fieldnames(pairs),'pairTransVectX' ), 1)))
    Angs=pairs(:,4:6);
else
    Angs=[pairs(:).pairTransAngleZXZPhi ; pairs(:).pairTransAngleZXZPsi ; pairs(:).pairTransAngleZXZTheta];
    Angs=Angs';
end;



meanAng=tom_average_rotations(Angs);
for i=1:size(Angs,1)
        lendiffAng(i)=tom_angular_distance(Angs(i,:),meanAng);
end;
%pdfAngOrg = fitdist(lendiffV','Kernel','BandWidth',4);
if (strcmp(angDist.prenorm,'1std'))
    lendiffAngN=lendiffAng./std(lendiffAng);
else
    lendiffAngN=lendiffAng;
end;

if (strcmp(angDist.type,'normal'))
    pdfA = makedist('normal','sigma',angDist.nStd.*std(lendiffAngN));
end;
if (vis.plotDistributions)
    x =0:0.1:max(lendiffAngN)*2.5;
    ySix = pdf(pdfA,x);
    figure;plot(x,ySix,'k-','LineWidth',2);
     smpTrain= pdf(pdfA,lendiffAngN);
    hold on; plot(lendiffAngN,smpTrain,'go'); hold off;
     title('pdf Ang');
end;

pdfST.AngAvg=meanAng;
pdfST.pdf=pdfA;
pdfST.data=lendiffAng;
pdfST.std=std(lendiffAng);
pdfST.prenorm=angDist.prenorm;

andDistMean=round(tom_angular_distance([0 0 0],meanAng(:)'));
meanAng=round(meanAng(:));

disp(['Angular dist mean: ' num2str(meanAng(1)) ' ' num2str(meanAng(2)) ' ' num2str(meanAng(3)) ' std: ' num2str(std(lendiffAng)) '   ang. Dist Mean: ' num2str(andDistMean) ]);
disp(' ');

function pdfST=genVectPdf(pairs,transDist,vis)

if ( isempty(find(ismember(fieldnames(pairs),'pairTransVectX' ), 1)))
    vects=pairs(:,1:3);
else
    vects=[pairs(:).pairTransVectX ; pairs(:).pairTransVectY ; pairs(:).pairTransVectZ];
    vects=vects';
end;

 meanV=mean(vects);
 diffV=vects-repmat(meanV,size(vects,1),1);
 for i=1:size(diffV,1)
      lendiffV(i)=norm(diffV(i,:));
end;
if (strcmp(transDist.prenorm,'1std'))
    lendiffVN=lendiffV./std(lendiffV);
else
    lendiffVN=lendiffV;
end;

if (strcmp(transDist.type,'normal'))
    pdfV = makedist('normal','sigma',transDist.nStd*std(lendiffVN));
end;

if (vis.plotDistributions)
    x =0:0.1:(max(lendiffVN))*2.5;
    ySix = pdf(pdfV,x);
    figure;plot(x,ySix,'k-','LineWidth',2);
    smpTrain= pdf(pdfV,lendiffVN);
    hold on; plot(lendiffVN,smpTrain,'go'); hold off;
    title('pdf Vect');
end;

pdfST.mean=meanV;
pdfST.pdf=pdfV;
pdfST.data=lendiffV;
pdfST.std=std(lendiffV);
pdfST.prenorm=transDist.prenorm;

meanV=round(meanV(:));
lenMenV=round(norm(meanV));
disp(['Vector   dist mean: ' num2str(meanV(1)) ' ' num2str(meanV(2)) ' ' num2str(meanV(3)) ' std: ' num2str(round(std(lendiffV)))  ' length meanVect: '  num2str(lenMenV) ]);


function pairsClean=cleanTrainData(pairs,transPairsParam,angPairsParam,vis)

disp(' ');
disp('  cleaning Distributions');


idxCleanTrans=cleanVectorDistribution(pairs,transPairsParam.clean,vis);
idxCleanAng=cleanAngleDistribution(pairs,angPairsParam.clean,vis);

idxCmb=intersect(find(idxCleanTrans),find(idxCleanAng));
pairsClean=pairs(idxCmb,:);

disp(' ');

function idxGood=cleanVectorDistribution(pairs,param,vis)
 
if ( isempty(find(ismember(fieldnames(pairs),'pairTransVectX' ), 1)))
    vects=pairs(:,1:3);
else
    vects=[pairs(:).pairTransVectX ; pairs(:).pairTransVectY ; pairs(:).pairTransVectZ];
    vects=vects';
end;

meanV=mean(vects(:,1:3));
idxGood=1:size(vects,1); 


if (vis.plotCleanUp==1)
    orig=zeros(size(vects));
    figure; quiver3(orig(:,2),orig(:,2),orig(:,3),vects(:,1),vects(:,2),vects(:,3),1,'color',[0 1 0]);    
    hold on; quiver3(0,0,0,meanV(:,1),meanV(:,2),meanV(:,3),1,'color',[0 0 1]);
    title('translation training data clean up');
end; 

for ii=1:param.iter
    meanV=mean(vects(idxGood,1:3));
    diffV=vects-repmat(meanV,size(vects,1),1);
    lendiffV=zeros(size(diffV,1),1);
    for i=1:size(diffV,1)
        lendiffV(i)=norm(diffV(i,:));
    end;
    thr=mean(lendiffV(idxGood))+(param.nStd*std(lendiffV));
    idxGood=lendiffV<thr;
    clear('lendiffV');
end;
disp(['translational cleaning removed: ' num2str(length(find(idxGood==0))) ' of '  num2str(length(idxGood)) ]);

if (vis.plotCleanUp==1)
    vectsBad=vects((idxGood==0),1:3);
    orig=zeros(size(vectsBad));
    quiver3(orig(:,2),orig(:,2),orig(:,3),vectsBad(:,1),vectsBad(:,2),vectsBad(:,3),1,'color',[1 0 0]);    
    quiver3(0,0,0,meanV(1),meanV(2),meanV(3),1,'color',[0 0 1]);
    plot3(meanV(1),meanV(2),meanV(3),'x','MarkerSize',25,'color',[0,0,0]);
    hold off;
end;


function pairs=alignToOneDirection(pairs)

vects=[pairs(:).pairTransVectX ; pairs(:).pairTransVectY ; pairs(:).pairTransVectZ];
vects=vects';
vectsInv=[pairs(:).pairInvTransVectX ; pairs(:).pairInvTransVectY ; pairs(:).pairInvTransVectZ];
vectsInv=vectsInv';


cl=kmeans(vects,2);

if (length(find(ismember(cl,1)))>length(find(ismember(cl,2))) )
    useCl=1;
else
    useCl=2;
end;
idx=find(cl==useCl);
meanV=mean(vects(idx,1:3));


diffV=vects-repmat(meanV,size(vects,1),1);
diffV=sqrt(sum(diffV(:,:).*diffV(:,:),2));
diffVInv=vectsInv-repmat(meanV,size(vectsInv,1),1);
diffVInv=sqrt(sum(diffVInv(:,:).*diffVInv(:,:),2));

idxSwap=find(diffV<diffVInv);

for i=1:length(idxSwap)
    pairs(idxSwap(i))=swapPairOrderEntry(pairs(idxSwap(i)));
end;

vects(idxSwap,:)=vectsInv(idxSwap,:);   
meanV=mean(vects);


function entry=swapPairOrderEntry(entry)


entry=swapFields(entry,'pairIDX1','pairIDX2');

%TransVect
entry=swapFields(entry,'pairTransVectX','pairInvTransVectX');
entry=swapFields(entry,'pairTransVectY','pairInvTransVectY');
entry=swapFields(entry,'pairTransVectZ','pairInvTransVectZ');

%TransAng
entry=swapFields(entry,'pairTransAngleZXZPhi','pairInvTransAngleZXZPhi');
entry=swapFields(entry,'pairTransAngleZXZPsi','pairInvTransAngleZXZPsi');
entry=swapFields(entry,'pairTransAngleZXZTheta','pairInvTransAngleZXZTheta');

%Coordinates
entry=swapFields(entry,'pairCoordinateX1','pairCoordinateX2');
entry=swapFields(entry,'pairCoordinateY1','pairCoordinateY2');
entry=swapFields(entry,'pairCoordinateZ1','pairCoordinateZ2');

%Angles
entry=swapFields(entry,'pairAnglePhi1','pairAnglePhi2');
entry=swapFields(entry,'pairAnglePsi1','pairAnglePsi2');
entry=swapFields(entry,'pairAngleTheta1','pairAngleTheta2');

%AddInfo
entry=swapFields(entry,'pairClass1','pairClass2');
entry=swapFields(entry,'pairPsf1','pairPsf2');





function st=swapFields(st,field1,field2)

v1=getfield(st,field1);
v2=getfield(st,field2);

st=setfield(st,field1,v2);
st=setfield(st,field2,v1);


function idxGood=cleanAngleDistribution(pairs,param,vis)


if ( isempty(find(ismember(fieldnames(pairs),'pairTransVectX' ), 1)))
    Angs=pairs(:,4:6);
else
    Angs=[pairs(:).pairTransAngleZXZPhi ; pairs(:).pairTransAngleZXZPsi ; pairs(:).pairTransAngleZXZTheta];
    Angs=Angs';
end;



AngAvg=tom_average_rotations(Angs);

if (vis.plotCleanUp)
    figure;
    tom_rotVectByAng([1 0 0],Angs,'zxz','vector',[0 1 0]);
    hold on;
    tom_rotVectByAng([1 0 0],AngAvg,'zxz','vector',[0 0 1]);
    title('rotation training data clean up');
end;

idxGood=1:size(Angs,1);
for ii=1:param.iter;
    AngAvg=tom_average_rotations(Angs(idxGood,:));
    for i=1:size(Angs,1)
       % lendiffV(i)=tom_angular_distance(AngAvg,Angs(idxGood(i),:));
        lendiffV(i)=tom_angular_distance(AngAvg,Angs(i,:));
    end;
    thr=mean(lendiffV(idxGood))+(param.nStd*std(lendiffV));
    idxGood=lendiffV<thr;
    disp([num2str(ii) ': thr: ' num2str(thr) ' nrPart: ' num2str(length(find(idxGood))) ' of ' num2str(size(Angs,1)) ]);
    clear('lendiffV');
end;
disp(['rotational cleaning removed: ' num2str(length(find(idxGood==0))) ' of '  num2str(size(Angs,1)) ]);


if (vis.plotCleanUp==1)
    AngsBad=Angs((idxGood==0),1:3);
    if (isempty(AngsBad)==0)
        tom_rotVectByAng([1 0 0],AngsBad,'zxz','vector',[1 0 0]);
    end;
    tmpAvg=tom_rotVectByAng([1 0 0],AngAvg,'zxz','vector',[0 0 1]);
    plot3(tmpAvg(1),tmpAvg(2),tmpAvg(3),'x','MarkerSize',25,'color',[0 0 0]);
    hold off;
end;




function transPairs=readTrainData(io,vis)

zz=1;
for i=1:length(io.trainList)
    [~,~,ext]=fileparts(io.trainList{1});
    if (strcmp(ext,'.star'))
           transPairs=tom_starread(io.trainList{i});
           if (io.useClass(1)>-1)
               headTmp=transPairs(1).Header;
               idx=find(ismember([transPairs(:).pairClass],io.useClass(i)));
               transPairs=transPairs(idx);
               transPairs(1).Header=headTmp;
           end;
    end;
        
    if (strcmp(ext,'.pair'))
        [~,~,ext]=fileparts(io.trainList{i});
        list=load(io.trainList{i});
        
        for ii=1:size(list,1)
            pair=list(ii,:);
            transPairs(zz,:)=calcPairTransForm(pair);
            zz=zz+1;
        end;
    end;

    if (vis.plotTrainInput)
            figure;
            tom_plot_vectorFied(io.trainList{i},[0 0 1],40,[0 0 1],1,'', io.useClass(i));
    end;
    
end;


