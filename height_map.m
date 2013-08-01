function height_map(c_lat,c_long,width)

[N,E]=GPS2UTM(c_lat,c_long);
N_b=N;E_b=E;
inc_lat=4e-4;
inc_acc=0;
while N_b<=N+width/2
    inc_acc=inc_acc+inc_lat;
    [N_b,E_b]=GPS2UTM(c_lat+inc_acc,c_long);
end
b_lat=c_lat+inc_acc;

inc_long=6e-4;
inc_acc=0;
while E_b<=E+width/2
    inc_acc=inc_acc+inc_long;
    [N_b,E_b]=GPS2UTM(c_lat,c_long+inc_acc);
end
b_long=c_long+inc_acc;

i=0;j=0;
for lat=2*c_lat-b_lat:inc_lat:b_lat
    i=i+1;j=0;
    for long=2*c_long-b_long:inc_long:b_long
        j=j+1;
        [N_temp,E_temp]=GPS2UTM(lat,long);
        land(i,j,1)=E_temp-E;
        land(i,j,2)=N_temp-N;
        land(i,j,3)=height_calculator(lat,long);
        pause(0.5);
        save('land','land','-mat'); 
    end
    i*j
end