% Function for calculate absolute and relative causal strengths  
% of a pair of time series
%
% Syntax: causal_matrix = causal_decomposition(s1,s2,rnoise,nensemble)
% Parameters:
% s1: time series 1
% s2: time series 2
% rnoise: the level of noise used in ensemble empirical mode decomposition
%         the rnoise value represents the fraction of standard deviation of
%         time series (e.g., 0.1)
% nensemble: the number of ensemble to average the results of noise-assisted
%            empirical mode decomposition. 
%            The error of IMFs equals to rnoise/sqrt(nensemble)        
% 
% Output:
% The output of causal_decomposition function contains a four by n matrix, 
% where n is the number of IMFs decomposed from the data. The first column 
% in the matrix indicates the relative causal strength from first time series 
% to the second time series and vise versa in the second column. The third 
% column indicates the absolute causal strength from the first time series 
% to the second time series, and vice versa in the fourth column. 
%
% Example: load ecosystem_data;
%          causal_matrix = causal_decomposition(DIDI,PARA,0.35,1000);
%
% Ver 1.0: Albert C. Yang, MD, PhD 7/11/2018
%
% Referece

function causal_matrix = causal_decomposition(s1,s2,rnoise,nensemble)

% normalize time series
s1=s1-mean(s1);
s1=s1/std(s1);
s2=s2-mean(s2);
s2=s2/std(s2);

% do EMD for each signal
if rnoise>0
    imfs1=eemd(s1,rnoise,nensemble,0)';
    imfs2=eemd(s2,rnoise,nensemble,0)';
else
    nensemble=1;
    imfs1=eemd(s1,rnoise,nensemble,0)';
    imfs2=eemd(s2,rnoise,nensemble,0)';
end

% determine imf size
imfsize=size(imfs1);
imfsize=imfsize(2)-1;

% initiate phase vectors
p12=zeros(imfsize,1);
p21=zeros(imfsize,1);

% calculate phase coherence between paired IMFs
[ps,pd]=phasefcimf(imfs1,imfs2);

% calculate variance of IMFs
v1=nvar(imfs1)';
v2=nvar(imfs2)';

% compute causal strength for each paired IMFs
for i=1:imfsize

   % remove an IMF and do the redecomposition
   s1r=s1-imfs1(:,i);
   s2r=s2-imfs2(:,i);
   imfs1r=eemd(s1r,rnoise,nensemble,0)';
   imfs2r=eemd(s2r,rnoise,nensemble,0)';  
   
   % recalculate phase coherence between paired IMFs
   ps12=phasefcimf(imfs1,imfs2r);
   ps21=phasefcimf(imfs1r,imfs2);
   
   % calculate absolute causal strength using variance-weighted Euclidian distance 
   % between the phase coherence of the original IMFs and redecomposed IMFs 
   p12(i)=sqrt(sum(v1.*v2.*(ps12-ps).^2)/sum(v1.*v2));
   p21(i)=sqrt(sum(v1.*v2.*(ps21-ps).^2)/sum(v1.*v2));
       
end

% calculate relative causal strengths from absolute causal strengths (p12 and p21)
% the output is 4 by n matrix (n: number of IMFs)
% the first column indicates relative causal strength from time series 1 to 2
% the second column indicates relative causal strength from time series 2 to 1
% the third column indicates absolute causal strength from time series 1 to 2
% the fourth column indicates absolute causal strength from time series 2 to 1
for i=1:imfsize
    
    if p12(i)<0.05 & p21(i)<0.05
        alpha=1;
    else
        alpha=0;
    end

    causal_matrix(i,:) = [(alpha+p12(i))/(2*alpha+p12(i)+p21(i)) (alpha+p21(i))/(2*alpha+p12(i)+p21(i)) p12(i) p21(i)];

end
