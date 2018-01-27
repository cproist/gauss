function [H, Pval, TestStat, CV, NullStat]= KSpr(Y,Nstraps)
% [H, Pval, TestStat, CV, NullStat]= prb(Y,Nphrands,Nboots)
% Phase-randomization bootstratp test for non-gaussianity/non-normality
% based on the Kolmogorov-Smirnov statistic.
%
% In: 
%  Y        = Time Series
%  Nstraps  = Number of phase randomized bootstraps
%
% Out:
%  H        = Hypothesis test - 1: Null Hypothesis of Gaussianity is rejected
%  Pval     = P-value
%  TestStat = Test Stastic
%  CV       = Critical Value
%  NullStat = Null Distribution. CV is computed as prctile(NullStat,95); 
%
%
% Reference :Proistosescu, C., P. Huybers, A. Rhines (2016), Identification
% and Interpretation of Non-Normality in Atmospheric Time-Series.
% Geophysical Research Letters, Vol 43 (Supplementary Information)
%
% cproist@gmail.com
% Last modified: 01/12/2018


% Ensure Input is a column 
if ~iscolumn(Y);
    Y=Y';
end;

% Create Nphhrands phase randomized iterations
[Ypr]=squeeze(phaseran(Y,Nstraps));    

% Compute Null-Hypothesis
NullStat=zeros(Nstraps,1);
mu=mean(Y);
sig=std(Y);
for ct=1:Nstraps
    y=(Ypr(:,ct)-mu)/sig;    
    binEdges    =  [-inf ; sort(y) ; inf];
    binCenters  =  (binEdges(1:end-1)+binEdges(2:end))*0.5;
    normCDF     =  normcdf(binCenters,0,1);            
    binCounts  =  histc (y , binEdges, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);    
    sampleCDF  =  sumCounts(1:end-1);        
    deltaCDF    =  abs(sampleCDF - normCDF);
    NullStat(ct)   =  max(deltaCDF);    
end
CV=prctile(NullStat,95);


% Compute Test Statistic
    y=(Y-mu)/sig;    
    binEdges    =  [-inf ; sort(y) ; inf];
    binCenters  =  (binEdges(1:end-1)+binEdges(2:end))*0.5;
    normCDF     =  normcdf(binCenters,0,1);
    binCounts  =  histc (y , binEdges, 1);
    sumCounts  =  cumsum(binCounts)./sum(binCounts);
    sampleCDF  =  sumCounts(1:end-1);
    deltaCDF    =  abs(sampleCDF - normCDF);
    TestStat   =  max(deltaCDF);    

TestStat=mean(TestStat);
Pval=sum((NullStat-TestStat)>0)/length(NullStat);
H=[TestStat>CV;Pval<0.05];

