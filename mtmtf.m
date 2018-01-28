% function [fq,H,ph] = mtmtf(x,y,dt,nw)
%
% Simple multi-taper method coherence estimate and Monte Carlo
% confidence estimation procedure. Based on Peter Huybers' multi-taper
% coherence estimator. http://www.people.fas.harvard.edu/~phuybers/Mfiles/
%
% In:
%      x       - invput vector 1.
%      y       - Input vector 2.
%      dt      - sampling interval
%      nw      - number of windows to use (8 is often a good choice)
%    
%
% Out:
%      fq      - frequency
%      H       - Transfer function (complex) 
%      ph      - phase
%
% Cristian Proistosescu
% Last modified 01/27/2018
% cproist@gmail.com


function [fq,H,ph] = mtmTF(x,y,dt,NW)
confn=0;
x   = x(:)-mean(x); y=y(:)-mean(y);

%check input
if nargin<6, qplot=0;  end;
if nargin<5, confn=0; end; 
if length(confn)==0, confn=0; end;
if nargin<4, NW=8;     end;
if length(NW)==0, NW=8;end;
if nargin<3, dt=1;     end;
if length(dt)==0, dt=1;end;

%define some parameters
N   = length(x);
k   = min(round(2*NW),N); 
k   = max(k-1,1);
fq   = (0:1/(N*dt):1/dt-1/(N*dt))';
%pls=2:(N+1)/2+1;
%if rem(length(y),2)==1; pls=pls(1:end-1); end;

%Compute the discrete prolate spheroidal sequences, requires the
%spectral analysis toolbox.
[E,V]=dpss(N,NW,k);

%Compute the windowed DFTs.
fkx=fft(E(:,1:k).*x(:,ones(1,k)),N);
fky=fft(E(:,1:k).*y(:,ones(1,k)),N);

%Compute coherence
Cxy= sum( (fkx.*conj(fky))' );
H  = sum((fkx.*conj(fky))')./sum([fkx.*conj(fkx)]');
ph = angle(Cxy)*180/pi;





