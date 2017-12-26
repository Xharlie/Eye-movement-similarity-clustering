function [BC, SC] = threshold(Cij)

% threshold.m - find proper threshold at which there is a phase transition
% for a network Cij; returns size, index, and mass of big and second
% clusters at steps of .99:0.01:-1.00

clear tF; clear triF;
tClusters = 1:length(Cij);
mass = ones(length(Cij),1)./length(Cij);
tF = triu(Cij,1);
[triF(:,1),triF(:,2)] = find(tF~=0);
triF(:,3) = tF((tF~=0));
sFij = flipdim(sortrows(triF,3),1); %Fij sorted in descending order
step = 1;
BC = [0 0 0];
SC = [0 0 0];

for thresh = 1.99:-0.01:-1
    incluster = find(sFij(:,3)>thresh);
    if incluster
        row = sFij(incluster,1);
        col = sFij(incluster,2);
        for j = 1:length(incluster)
            [parentrow,parentcol] = find(tClusters==row(j));
            [addrow,addcol] = find(tClusters==col(j));
            if addcol~=parentcol
                addons = tClusters(find(~ismember(tClusters(:,addcol),tClusters(:,parentcol))),addcol);
                tClusters((end+1):(end+length(addons)),parentcol) = addons;
                tClusters(:,addcol) = 0;
                mass(parentcol) = mass(parentcol)+mass(addcol);
                mass(addcol) = 0;
            end
        end
    end
    for c=1:length(Cij)
        cLength(c) = length(find(tClusters(:,c)~=0));
    end
    if cLength(find(cLength>1))
        BC(step,:) = [max(cLength(find(cLength>1))) min(find(cLength==max(cLength(find(cLength>1))))) mass(min(find(cLength==max(cLength(find(cLength>1))))))];
        nextClusters = cLength(cLength~=max(cLength(find(cLength>1)))); %find next-largest clusters
        if nextClusters(find(nextClusters>1))
            SC(step,:) = [max(nextClusters(find(nextClusters>1))) min(find(nextClusters==max(nextClusters(find(nextClusters>1))))) mass(min(find(nextClusters==max(nextClusters(find(nextClusters>1))))))];
        end
    end
    step = step+1;
end
if length(SC(:,1))<length(BC(:,1))
    SC(length(BC(:,1)),3) = 0;
end

% plot(1.99:-0.01:-1,BC(:,1)); hold on; plot(1.99:-0.01:-1,SC(:,1),'r'); hold off;
% title(['Networks built up according to correlations']);
% xlabel('Threshold'); ylabel('Mass of cluster'); legend('Largest cluster','Second-largest cluster');
% axis([-1 2 0 27]);
% (200-max(find(diff(BC(:,1))>0))-1)/100
% pause
%         grtitle = cellstr(['C:/Users/Ka/Documents/JAMLab/F_ij/clustered_Fij/clusterEvolution', num2str(vid),'.png']);
%         saveas(figure(1),grtitle{1});
%         clf;
end
