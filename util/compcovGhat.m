% Function which determines the estimate of cov(Ghat) in the frequency domain (see slide 2-24)
% This function only works when the model has been identified in a model structure M containing both Go and Ho
%
% Use:
%
% [covG,G]=compcovGhat(sys,w)
%
% Input arguments: 
%
% "sys" is the identified model expressed in the IDPOLY-object format
% "w" is a vector containing the frequencies at which you want cov(hatG) to be evaluated
%
% Output arguments:
%
% "covG" is a vector containing cov(Ghat) evaluated at the frequencies in "w"
% "G" is a vector containing the modulus of the identified model Ghat evaluated at the frequencies in "w"

function [covG,G]=compcovGhat(sys,w)




[Gi,p,gccov]=freqresp(sys,w);
nwdef=length(w);

vect2=[]; G=[];
for ii=1:nwdef,
    mat=squeeze(gccov(:,:,ii,:,:));
vect2=[vect2,(trace(mat))]; G=[G,Gi(:,:,ii)];
end;
G=abs(G);
covG=vect2;



% [mag,ph,ww,st]=boderesp(sys,w);
% 
% vect=[]; 
% for ii=1:nwdef,
%     
% vect=[vect,st(:,:,ii)];
% end;
