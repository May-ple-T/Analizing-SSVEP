% Demo of MsetCCA for four-class SSVEP Recognition in BCI %
% by Yu Zhang, ECUST, 2013.6.12
% Email: zhangyu0112@gmail.com

clc;
clear all;
close all;


%% Initialize parameters
Fs=128;                                  % sampling rate
t_length=5;                              % data length (4 s)
TW=1:1:t_length;
TW_p=round(TW*Fs);
n_run=20;                                % number of used runs
sti_f=[15 12 10 8.6];             % stimulus frequencies 15, 12, 10, 8.57 Hz
n_sti=length(sti_f);                     % number of stimulus frequencies
n_correct=zeros(3,length(TW));


%% Load SSVEP data
load obabw.mat
% Data description:
% 8 channels x 512 points x 20 trials x 4 stimulus frequencies


%% CCA for SSVEP recognition
% Construct reference signals of sine-cosine waves
N=2;    % number of harmonics
ref1=refsig(sti_f(1),Fs,t_length*Fs,N);
ref2=refsig(sti_f(2),Fs,t_length*Fs,N);
ref3=refsig(sti_f(3),Fs,t_length*Fs,N);
ref4=refsig(sti_f(4),Fs,t_length*Fs,N);

% Recognition
for run=1:n_run
    for tw_length=1:t_length       % time window length:  1s:1s:4s
        fprintf('CCA Processing... TW %fs, No.crossvalidation %d \n',TW(tw_length),run);
        for j=1:n_sti
            [wx1,wy1,r1]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),ref1(:,1:TW_p(tw_length)));
            [wx2,wy2,r2]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),ref2(:,1:TW_p(tw_length)));
            [wx3,wy3,r3]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),ref3(:,1:TW_p(tw_length)));
            [wx4,wy4,r4]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),ref4(:,1:TW_p(tw_length)));
            [v,idx]=max([max(r1),max(r2),max(r3),max(r4)]);
            if idx==j
                n_correct(1,tw_length)=n_correct(1,tw_length)+1;
            else
                ans=j;
            end
        end
    end
end


%% MsetCCA for SSVEP recognition
K=2;    % number of extracted components for each spatial filter
for run=1:n_run
    idx_traindata=1:n_run;
    idx_traindata(run)=[];
    for tw_length=1:t_length       % time window length:  1s:1s:4s
        fprintf('MsetCCA Processing... TW %fs, No.crossvalidation %d \n',TW(tw_length),run);
        % Reference signals optimization by MsetCCA
        Temp1=zeros((n_run-1)*K,TW_p(tw_length)); Temp2=Temp1; Temp3=Temp2;  Temp4=Temp3;
        W1=msetcca(SSVEPdata(:,1:TW_p(tw_length),idx_traindata,1),K);
        W2=msetcca(SSVEPdata(:,1:TW_p(tw_length),idx_traindata,2),K);
        W3=msetcca(SSVEPdata(:,1:TW_p(tw_length),idx_traindata,3),K);
        W4=msetcca(SSVEPdata(:,1:TW_p(tw_length),idx_traindata,4),K);
        for qq=1:n_run-1
            Temp1((qq-1)*K+1:qq*K,:)=W1(:,:,qq)'*SSVEPdata(:,1:TW_p(tw_length),idx_traindata(qq),1);
            Temp2((qq-1)*K+1:qq*K,:)=W2(:,:,qq)'*SSVEPdata(:,1:TW_p(tw_length),idx_traindata(qq),2);
            Temp3((qq-1)*K+1:qq*K,:)=W3(:,:,qq)'*SSVEPdata(:,1:TW_p(tw_length),idx_traindata(qq),3);
            Temp4((qq-1)*K+1:qq*K,:)=W4(:,:,qq)'*SSVEPdata(:,1:TW_p(tw_length),idx_traindata(qq),4);
        end
        
        % Recognition
        for j=1:n_sti
            [wx1,wy1,r1]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Temp1);
            [wx2,wy2,r2]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Temp2);
            [wx3,wy3,r3]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Temp3);
            [wx4,wy4,r4]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Temp4);
            [v,idx]=max([max(r1),max(r2),max(r3),max(r4)])
            if idx==j
                n_correct(2,tw_length)=n_correct(2,tw_length)+1;
            end
        end
    end
end


%% IT-CCA for SSVEP recognition
% Construct reference signals
Chi=it_cca(SSVEPdata);
for j=1:n_sti
    [wx1,wy1,r1]=cca(Chi(:,1:TW_p(tw_length),j),ref1(:,1:TW_p(tw_length)));
    [wx2,wy2,r2]=cca(Chi(:,1:TW_p(tw_length),j),ref2(:,1:TW_p(tw_length)));
    [wx3,wy3,r3]=cca(Chi(:,1:TW_p(tw_length),j),ref3(:,1:TW_p(tw_length)));
    [wx4,wy4,r4]=cca(Chi(:,1:TW_p(tw_length),j),ref4(:,1:TW_p(tw_length)));
    [v,idx]=max([max(r1),max(r2),max(r3),max(r4)])
end

% Recognition
for run=1:n_run
    for tw_length=1:t_length       % time window length:  1s:1s:4s
        fprintf('IT-CCA Processing... TW %fs, No.crossvalidation %d \n',TW(tw_length),run);
        for j=1:n_sti
            [wx1,wy1,r1]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Chi(:,1:TW_p(tw_length),1));
            [wx2,wy2,r2]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Chi(:,1:TW_p(tw_length),2));
            [wx3,wy3,r3]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Chi(:,1:TW_p(tw_length),3));
            [wx4,wy4,r4]=cca(SSVEPdata(:,1:TW_p(tw_length),run,j),Chi(:,1:TW_p(tw_length),4));
            [v,idx]=max([max(r1),max(r2),max(r3),max(r4)]);
            [vmin,idxmin]=min([max(r1),max(r2),max(r3),max(r4)]);
            if idx==j
                n_correct(3,tw_length)=n_correct(3,tw_length)+1;
            end
        end
    end
end


%% Plot accuracy
accuracy=100*n_correct/n_sti/n_run
col={'b-*','r-o','g-x'};
for mth=1:3
    plot(TW,accuracy(mth,:),col{mth},'LineWidth',1);
    hold on;
end
xlabel('Time window length (s)');
ylabel('Accuracy (%)');
grid;
xlim([0.75 4.25]);
ylim([0 100]);
set(gca,'xtick',1:5,'xticklabel',1:5);
% title('\bf IT-CCA vs CCA for SSVEP Recognition');
% h=legend({'CCA','IT-CCA'});
title('\bf CCA for SSVEP Recognition');
h=legend({'CCA'});
set(h,'Location','SouthEast');

