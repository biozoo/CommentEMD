% Causal decompostion of moran effect time series
% Causal decomposition and eemd functions are from Yang et al. 2018 Nature Commun. and Wang et al. 2014 PhysicaA, respectively)
% load Moran time series data produced in R(Please change the route for dataset and functions before run this code). 
% Please add the EMD functions deposited in the file, 'causal_decomposition_20180901', to the path (selected folders & subfolders) 
load 'N1.txt' 
load 'N2.txt'   
load 'time_moran.txt'   

N1s=(N1-mean(N1))/std(N1);
N2s=(N2-mean(N2))/std(N2);
n=size(N1s,1)
plot_pairedimfs(time_moran,N1s,N2s,0.35,200,[]);

% preliminary causal decomposition
causal_matrix = causal_decomposition(N1s,N2s,0.35,200)
nimf=floor(log2(n))


% OI and separability test for N1
NE = 1000; % # of ensemble
numImf = nimf; % # of imfs
data=N1s
testr = 0.01:0.01:1
nr=size(testr,2)
rms=zeros(1,nr)
nonorth=zeros(1,nr)
for i = 1:nr
	Nstd = testr(i)
	[imf] = eemd(data,Nstd,NE,numImf)
	[rho, pval] = corr(transpose(imf), 'rows','pairwise');
	%%RMS in separability test
	U = triu(rho,1)
	rhoall = U(U~=0)
	rms(i)=mean(rhoall.^2)^0.5
	%%OI in Orthogonality test
	nonorth(i)=io(data,imf)	
	
end
emsem_r1=testr(rms==min(rms)==1)

% OI and separability test for N2
numImf = nimf; % # of imfs
data=N2s
testr = 0.01:0.01:1
nr=size(testr,2)
rms=zeros(1,nr)
nonorth=zeros(1,nr)
for i = 1:nr
	Nstd = testr(i)
	[imf] = eemd(data,Nstd,NE,numImf)
	[rho, pval] = corr(transpose(imf), 'rows','pairwise');
	%%RMS in separability test
	U = triu(rho,1)
	rhoall = U(U~=0)
	rms(i)=mean(rhoall.^2)^0.5
	%%OI in Orthogonality test
	nonorth(i)=io(data,imf)	
end

emsem_r2=testr(rms==min(rms)==1)



% Fig. 1a run causal_decomposition using updated emsemble parameters  
emsem_r=(emsem_r1+emsem_r2)/2
causal_matrix = causal_decomposition(N1s,N2s,emsem_r,1000)
figure
h=bar(causal_matrix(:,1:2),'stacked','b','r')
set(h,{'FaceColor'},{'b';'r'});
hline = refline(0,0.5);
hline.Color = 'k';
hline.LineStyle = '--';


