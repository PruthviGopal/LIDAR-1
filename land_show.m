function [centr_alt,minz,maxz]=land_show

close all;
h=load('land','-mat');

photo=zeros(size(h.land(:,:,:)));

maxz=max(max(h.land(:,:,3)));minz=min(min(h.land(:,:,3)));
photo(:,:,1)=( (h.land(:,:,3)-minz) / (maxz-minz) );
photo(:,:,2)=( (h.land(:,:,3)-minz) / (maxz-minz) );
photo(:,:,3)=( (h.land(:,:,3)-minz) / (maxz-minz) );

surf(h.land(:,:,1),h.land(:,:,2),h.land(:,:,3)...
    ,photo);
shading interp

centr_alt=h.land(ceil(size(h.land,2)/2),ceil(size(h.land,3)/2),3);