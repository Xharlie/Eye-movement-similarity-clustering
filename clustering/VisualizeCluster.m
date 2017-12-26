clear;
addpath('thermal/')
data_dir='../result/30s/c_location/';
data_dir_pattern='../result/30s/c_location/*lam*';
result_dir ='clustering_result/30s/c_location/'
Subs_list_c = [2,4,5,6,7,9,10,11,12,13,14,15,16;
               2,4,5,6,7,0,10,11,12,13,0,0,16;
               2,4,5,6,7,8,10,11,12,13,0,15,16;
               2,4,5,6,7,0,10,11,12,13,0,15,16;
               2,4,5,6,7,8,10,11,12,0,0,0,16;
               2,4,5,6,7,0,10,11,12,0,0,0,16;
               2,4,5,6,7,8,10,11,12,0,0,0,16;
               2,3,4,5,6,7,8,10,11,12,14,15,16;
               2,4,5,6,7,9,10,11,12,0,0,0,0;
               2,4,5,6,7,9,10,11,12,13,0,0,16;
               2,4,5,6,0,0,10,11,12,0,14,0,16;
               2,4,5,6,7,0,10,11,12,0,14,0,16;
               2,4,5,6,0,0,10,11,12,13,14,15,16;
               2,4,5,6,7,0,10,11,12,13,0,15,16;
               2,4,0,6,0,0,10,11,12,0,0,15,16;
               2,4,5,6,7,0,10,11,12,0,0,0,16;
               2,3,4,5,6,0,0,10,11,12,13,0,16;
               2,3,4,5,6,7,0,10,11,12,0,14,16;
               2,3,4,5,6,7,8,10,11,12,13,0,16;
               2,3,4,5,6,7,0,10,11,12,13,15,16;];
Label=["an1", "at1", "azz1", "cd1", "clt", "lr1", "lv1", ...
    "mcy1", "msb", "pc1", "pp1", "rn1", "sa", "sd2", "sw1", "ys1"];
n=size(Subs_list_c,1)
r=0;
NumSubsTotal=max(max(Subs_list_c));
Cluster=zeros(n,NumSubsTotal);
%Cluster = load([Dir, '/c_Cluster_r', num2str(r), '.mat'])
sfile_handle = dir('./Data/clean_csv_files_with_NANs/*_c_cleaner.csv');
data_sub_dir = dir(data_dir_pattern);
dirFlags = [data_sub_dir.isdir];
data_sub_dir=data_sub_dir(dirFlags);
for di = 1:size(data_sub_dir)
    for videoID=1:n
        file=strcat(data_sub_dir(di).name,'/dis_',int2str(videoID),'.mat')
        mkdir(strcat(result_dir,data_sub_dir(di).name));
        Node={}
        Edge={}
        Subs = Subs_list_c(videoID,:);
        Subs = Subs(Subs ~=0 )
        %Jinf = textread([Dir, '/Jinf_v0_video',num2str(videoID),'.txt']);
        %Jinf = Jinf(1:length(Jinf)-1,2);
    %     Jinf = textread([Dir, '/corr_corr_v0_video',num2str(videoID),'.txt']);
    %     Jinf = Jinf(1:length(Jinf),2);
        Jij = load(strcat(data_dir,file));
        Jij=Jij.matrix;
        Jij = squeeze(Jij)
        k=1;


    %     for i=1:length(Subs)-1
    %         for j=i+1:length(Subs)
    %             Jij(i,j) = Jinf(k);
    %             k=k+1;
    %         end
    %     end
    %     Jij = Jij+Jij';
    % 
    %     for i=1:length(Subs)
    %         Jij(i,i)=1;
    %     end
        [BC, SC] = threshold(Jij);
        %prompt='Threshold?'
        %PercolationThreshold=input(prompt);
        PercolationThreshold = (200-max(find(diff(BC(:,1))>0))-1)/100
        Jij(find(abs(Jij)<PercolationThreshold))=0;

        [coms,Qs] = mscd_afg(Jij, r);
        com=coms{1,1}
        for i=1:length(Subs)
            Cluster(videoID,Subs(i))=com(i);
        end

        Node{1,1}='Id';
        Node{1,2}='Label';
        Node{1,3}='Category'
        for i=1:length(Subs)
            Node{i+1,1}=Subs(i);
            Node{i+1,2}=Label(Subs(i))
            Node{i+1,3}=Cluster(videoID, Subs(i))
        end

        Edge{1,1}='Source';
        Edge{1,2}='Target';
        Edge{1,3}='Type';
        Edge{1,4}='Weight';
        %ClusterTemp=Cluster.Cluster(videoID,:);
        EdgeIndex=2;

        for i=1:(length(Jij)-1)
            for j=(i+1):length(Jij)
                if Jij(i,j)>0
                    Edge{EdgeIndex,1}=Subs(i);
                    Edge{EdgeIndex,2}=Subs(j);
                    Edge{EdgeIndex,3}='Undirected';
                    Edge{EdgeIndex,4}=Jij(i,j)
                    EdgeIndex=EdgeIndex+1;
                end
            end
        end


    %     for i=1:(length(ClusterTemp)-1)
    %         if ClusterTemp(i)>0
    %             list=find(ClusterTemp==ClusterTemp(i))
    %             for j=1:length(list)
    %                 if list(j)>i
    %                     Edge{EdgeIndex,1}=i;
    %                     Edge{EdgeIndex,2}=list(j);
    %                     Edge{EdgeIndex,3}='Undirected';
    %                     Edge{EdgeIndex,4}=Jij(find(Subs==i),find(Subs==list(j)))
    %                     EdgeIndex=EdgeIndex+1;
    %                 end
    %             end
    %         end
    %     end
        cell2csv([result_dir, data_sub_dir(di).name,'/corr_Node_video', num2str(videoID), '.csv'], Node);
        cell2csv([result_dir, data_sub_dir(di).name,'/corr_Edge_video', num2str(videoID), '.csv'], Edge);
    end
end