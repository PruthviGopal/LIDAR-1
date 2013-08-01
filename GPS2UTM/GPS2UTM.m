function [N,E]=GPS2UTM(lat,long)

[N,E,zone,lcm]=ell2utm(deg2rad(lat),deg2rad(long));