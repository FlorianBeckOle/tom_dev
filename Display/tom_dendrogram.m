function [groups,cmap,groupIdx,ColorThreshold]=tom_dendrogram(tree,ColorThreshold,nrObservations,dsp,maxLeaves)
%TOM_DENDROGRAM creates a dendrogram from linked data and gives the members
%and color for clusters
%
%   [groups,cmap,groupIdx]=tom_dendrogram(tree,ColorThreshold,nrObservations,dsp)
%
%PARAMETERS
%
%  INPUT
%   tree                          linked data
%   ColorThreshold         (mean(tree(:,3))) threshold for creating clusters
%   nrObservations         (-1) number of obs b4 linkage needed for cluster
%                                           members
%   dsp                            (1) display flag
%   maxLeaves               (2500) max number of leaves in dendrogram use 0
%                                                  to switch off
%   
%
%  OUTPUT
%   groups           		     structure containing members,color and id for
%                                       each cluster
%   cmap                          used colors
%   groupIdx                    index of groups
%   ColorThreshold         threshold applyed
%
%EXAMPLE
%
%
%
%REFERENCES
%
%SEE ALSO
%   linkage,dendrogram
%
%
%   updated by FB 09/08/19
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
    ColorThreshold=max(tree(:,3))*0.7;
end;

if (nargin<3)
    nrObservations=-1;
end;

if (nargin<4)
    dsp=1;
end;

if (nargin<5)
    maxLeaves=1500;
end;

if (ischar(tree))
    tree=tom_mrcread(tree);
    tree=tree.Value;
end;




if (strcmp(ColorThreshold,'auto'))
    clear(ColorThreshold)  ;
    ColorThreshold=max(tree(:,3))*0.7;
end;

if (nrObservations>maxLeaves)
    if (dsp==1)
        set(gcf,'Name',['clustering '  num2str(maxLeaves) ' of ' num2str(nrObservations) ' shown']);
    end;
end;


[groupIdx,~,cmap]=genColorLookUp(tree,ColorThreshold);

groups='';
dlabels='';
ugroupIdx=unique(groupIdx);
for i=1:length(ugroupIdx)
    groups(i).id=ugroupIdx(i);
    if (groups(i).id==0)
        groups(i).color=cmap(end,:);
    else
        groups(i).color=cmap(groups(i).id,:);
    end;
   tt=find(groupIdx==groups(i).id);
   tmpMem=unique(tree(tt,1:2));
   if (nrObservations>0)
        tmpMem=tmpMem(find(tmpMem<= nrObservations));
        if (isempty(tmpMem))
            tmpMem=-1;
        end;
   end;
   groups(i).members=tmpMem;
   groups(i).tree=tree(tt,:);
   if (groups(i).members>-1)
       for iii=1:length(groups(i).members)
           dlabels{groups(i).members(iii)}=['c' num2str(groups(i).id)];
       end;
   end;
end;

if (dsp)
  
    %label{i+1}='c0';
    if (isempty(dlabels))
        dendrogram(tree,maxLeaves,'ColorThreshold',ColorThreshold);
    else
        dendrogram(tree,maxLeaves,'ColorThreshold',ColorThreshold,'Labels',dlabels);
    end;
    %dendrogram(tree,0,'ColorThreshold',ColorThreshold);
end;


if (dsp && isempty(groups)==0)
    hold on;
    for i=1:length(groups)
        h(i)=plot(1,1,'color',groups(i).color);
        lText{i}=['cl:' num2str(groups(i).id) '(' num2str(length(groups(i).members)) ')'];
    end;
    legend(h,lText,'Location','westoutside','FontSize',10);
    hold off;
end;

function [theGroups,groups,cmap]=genColorLookUp(Z,threshold)

    if (strcmp(threshold,'off'))
        theGroups='';
        groups='';
        cmap='';
    else
        Z = transz(Z);
        
        theGroups = 1;
        groups = 0;
        cmap = [1 0 0];
        numLeaves = size(Z,1)+1;
     
        
        groups = sum(Z(:,3)< threshold);
        if groups > 1 && groups < (numLeaves-1)
            theGroups = zeros(numLeaves-1,1);
            numColors = 0;
            for count = groups:-1:1
                if (theGroups(count) == 0)
                    P = zeros(numLeaves-1,1);
                    P(count) = 1;
                    P = colorcluster(Z,P,Z(count,1),count);
                    P = colorcluster(Z,P,Z(count,2),count);
                    numColors = numColors + 1;
                    theGroups(logical(P)) = numColors;
                end
            end
            cmap = hsv(numColors);
            cmap(end+1,:) = [0.7 0.7 0.7];
        else
            theGroups='';
            groups='';
            cmap='';
        end
    end;
    
   
 function T = colorcluster(X, T, k, m)
% find local clustering

n = m;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = colorcluster(X, T, k, n); % trace back left subtree
        T = colorcluster(X, T, X(n,2), n);
        break;
    end
end
T(m) = 1;

function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

numLeaves = size(Z,1)+1;

for i = 1:(numLeaves-1)
    if Z(i,1) > numLeaves
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > numLeaves
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end

function a = traceback(Z,b)

numLeaves = size(Z,1)+1;

if Z(b-numLeaves,1) > numLeaves
    a = traceback(Z,Z(b-numLeaves,1));
else
    a = Z(b-numLeaves,1);
end
if Z(b-numLeaves,2) > numLeaves
    c = traceback(Z,Z(b-numLeaves,2));
else
    c = Z(b-numLeaves,2);
end

a = min(a,c);

   