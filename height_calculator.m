function h=height_calculator(lat,long)

to_google=['http://maps.googleapis.com/maps/api/elevation/json?locations='...
    num2str(lat) ',' num2str(long) '&sensor=false'];
a=urlread(to_google);
h=str2num(a(regexp(a,'"elevation" : ','end')+1:regexp(a,',')-1));