%for long videos that need to be split up
% generate required format if not already in memory or disk
if ~exist('X')
    if ~exist('data_myformat_split.mat') 
        
        load full_eye_data %just change the path to wherever you have 'data.m' stored
        
        s=data(:,9); % subject
        v=data(:,10); % video
        x=data(:,3); % horixontal eye position
        y=data(:,5); % vertical eye position
        a=data(:,7); % pupil area
        S = max(s);
        V = [1 3 5 6 9:13 18 20:22 27 29 32 36 37 40:42 45 46 51 56 59 60];
        fs=250;
        idx = 1;
        
        clear X;
        for n=1:length(V)
            i = V(n);
            for j=1:max(s)
                indx = find(v==i&s==j);
                if ~isempty(indx)
                    t = 0;
                    if ismember(i, [1 3 5 6 9:13 18 20 21 27 29 32 36 37 40:42 51 56 59 60])
                        t(2)=length(indx)-7250;
                        t(3)=length(indx);
                    end
                    
                    if ismember(i, [22 46])
                        t(2)=length(indx)-14500;
                        t(3)=length(indx)-t(2);
                        t(4)=length(indx);
                    end
                    
                    if i==45
                        t(2)=length(indx)-21750;
                        t(3)=length(indx)-(2*t(2));
                        t(4)=length(indx)-t(2);
                        t(5)=length(indx);
                    end
                    for m=2:length(t)
                        X{idx+m-2}(:,j) = x(indx(t(m-1)+1:t(m)));
                        Y{idx+m-2}(:,j) = y(indx(t(m-1)+1:t(m)));
                        A{idx+m-2}(:,j) = a(indx(t(m-1)+1:t(m)));
                        xqual(idx+m-2,j) = mean(data(indx(t(m-1)+1:t(m)),4));
                        vidList(idx+m-2,1:2) = [i idx+m-1];
                    end
                end
            end
            % show us a little something to get an idea
            subplot(3,1,1); plot(A{idx});
            title(['pupil area, Video #' num2str(idx)]);
            subplot(3,1,2); plot(X{idx});
            title(['vertical eye position, Video #' num2str(idx)]);
            subplot(3,1,3); plot(Y{idx});
            title(['horizontal eye position, Video #' num2str(idx)]);
            drawnow
            idx = length(vidList(:,2))+1;
        end
        save data_myformat_split.mat X Y A S V xqual fs;
    else
        load data_myformat_split.mat
    end
end

% these are subjects with apparent delay (form visual inspection)
delay = [];
delay(:,1)=[2 12 24 37 15 18 33 35 29]'; % these subjects seem delayed.
delay(:,2)=[1  1  1  1  1  1  1  1 5]'; % by this much

% subjects we have to remove (even if they slip though automatic removal)

load subs_remove; % I went through all of the plots a few times and picked out the bad viewers, then made a list
%subs_remove{i} = []; 

% window to remove samples before and after nan (blinks) (20ms)
mask = ones(round(0.04*fs),1);

%A{46}(:,18)=nan; % not sure how else to get rid of this one

feature={'Pupil area','Horizontal eye position','Vertical eye position'};
F=length(feature);

addpath('/home/lingchen/Desktop/Kates_Stuff/MatlabCode/inexact_alm_rpca/')
addpath('/home/lingchen/Desktop/Kates_Stuff/MatlabCode/inexact_alm_rpca/PROPACK') %adding paths of code to clean data - change if necessary

for i=1:size(vidList,1) % videos
    
    for j=1:F % features
        
        switch j
            case 1, x=A{i}; crange=[0 1000];
            case 2, x=X{i}; crange=[0 2400];
            case 3, x=Y{i}; crange=[0 2400];
        end
        
        % shift by delay
        for k=1:size(delay,1)
            d=round(delay(k,2)*fs); s=delay(k,1);
            x(:,s)=[nan(d,1); x(1:end-d,s)];
        end
        
        % keep subjects with more than certain porcentate of good data
        subs = find(mean(isnan(x)|x==0)<1-0.8);
        
        % remove specific problem subjects
        [~,indx,~] = intersect(subs,subs_remove{vidList(i,1)}); subs(indx)=[]; % KATE change {1}->{i}
        
        % cut out maximum delay, keep only subs, zero nan
        d = max(round(delay(:,2)*fs));
        x=x(d+1:end,subs);
        
        % axis labels
        [T,N] = size(x); T=(1:T)/fs; N=1:N;
        Tkeeper_split{i,j} = T;
        Nkeeper_split{i,j} = N;
        
        % show the data
        subplot(F,2,2*j-1); imagesc(T,N,x')
        if j==F, xlabel('Time (s)'); end
        title([ feature{j} ', Video #' num2str(vidList(i,1))]);
        ylabel('subjects');
        caxis(crange)
        
        % remove samles flanking nan
        x=flipud(filter(mask,1,flipud(filter(mask,1,x)))); 
        x=x(length(mask):end-length(mask)+1,:); % remove edge we messed up
        
        %and now clean the data using a sparse PCA method
        lambda = 0.015;

        x(isnan(x))=0; 
        [x_clean x_outlier] = inexact_alm_rpca(x,lambda);
        subplot(F,2,2*j); imagesc(T,N,x_clean');
        title('cleaned')
        caxis([0 2400])
        saveR_split{i,j} = x;
        
        makeR_split{i,j}=x_clean; %this is where the cleaned data is stored; the 
        %x and y data is stored separately so you'll need to make 'r' out
        %of sqrt(x^2 + y^2) and find the correlation coefficient after the
        %whole loop completes
        
        % show correlation coefficient
        R=nancorrcoef(x_clean);
        subplot(F,2,j*2); imagesc(R); axis equal; axis tight; colorbar
        if j==1, title('Correlation coefficients'); end 
        
        cc(i,j) = mean(mean(R-eye(max(N))));
        
        if isnan(cc(i,j)), keyboard; end
        
    end % features
    
%     drawnow;
%     grtitle = cellstr(['/home/kate/EyeData/cleanedPlots/vid', num2str(i),'.png']);
%     saveas(figure(1),grtitle{1});
    %pause
    
end % videos


figure(2); clf
load master2.mat % contains ratings and video names
for i=1:length(master2)
    if ~isempty(master2(i).Rating)
        ratings(i,1)=master2(i).Rating;
    else
        ratings(i,1)=nan;
    end
end

indx=find(~isnan(ratings))';
[c,p]=corrcoef([ratings(indx) cc(indx,:)]);
c=c(2:end);p=p(2:end);

for j=1:F
  subplot(2,F,j); plot(cc(indx,j),ratings(indx),'*');
  title(['r=' num2str(c(j),2) ', p=' num2str(p(j),2) ', N=' num2str(length(indx))]);
  xlabel(feature{j})
  ylabel('AdMeter rating')
  axis square
end

load ratingsSubjects65ads %there are two different sets of individual ratings; 
%this is the one with fewer ratings 
ratings=nanmean(ratingsSubjects65ads,2);
[c,p]=corrcoef([ratings cc]);
c=c(2:end);p=p(2:end);

for j=1:F
  subplot(2,F,F+j); plot(cc(:,j),ratings,'*');
  title(['r=' num2str(c(j),2) ', p=' num2str(p(j),2) ', N=' num2str(length(cc))]);
  xlabel(feature{j})
  ylabel('AdMeter rating')
  axis square
end

saveas(2,'correlations.jpg')














return

addpath('/home/kate/MatlabCode/inexact_alm_rpca')
addpath('/home/kate/MatlabCode/inexact_alm_rpca/PROPACK')

% left-over code to clean using sparce PCA.

%Donald - I put this in the body of the code earlier; you don't need to
%worry about what's below, it was just for earlier when the code was being
%written
lambda = 1 / sqrt(length(subs));
lambda = 0.015;


a(isnan(a))=0;
[a_clean a_outlier] = inexact_alm_rpca(a,lambda);
subplot(3,2,2); imagesc(T,N,a_clean')
title('cleaned')
caxis([0 1000])

x(isnan(x))=0;
[x_clean x_outlier] = inexact_alm_rpca(x,lambda);
subplot(3,2,4); imagesc(T,N,x_clean');
title('cleaned')
caxis([0 2400])

y(isnan(y))=0;
[y_clean y_outlier] = inexact_alm_rpca(y,lambda);
subplot(3,2,6); imagesc(T,N,y_clean')
title('cleaned')
caxis([0 1600])
drawnow








