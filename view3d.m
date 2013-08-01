function view3d(file_name,gate_size,SNR_lim,vel_lim,centr_alt,minz,maxz)

a=vel_load(file_name);

for beam=1:size(a,1)/(max(a(:,1))+1)
    for gate=1:max(a(:,1))+1
        row=(beam-1)*(max(a(:,1))+1)+gate;
        if a(row,3)>SNR_lim & a(row,2)<vel_lim %& (a(row,6)==17.2)
            x(gate,beam)=gate_size*a(row,1)*cosd(a(row,6))*sind(a(row,5));
            y(gate,beam)=gate_size*a(row,1)*cosd(a(row,6))*cosd(a(row,5));
            z(gate,beam)=gate_size*a(row,1)*sind(a(row,6))+centr_alt;
            c(gate,beam)=a(row,2);
        else
            x(gate,beam)=nan;
            y(gate,beam)=nan;
            z(gate,beam)=nan;
            c(gate,beam)=nan;
        end
        if z(gate,beam)>200
            x(gate,beam)=nan;
            y(gate,beam)=nan;
            z(gate,beam)=nan;
            c(gate,beam)=nan;
        end
    end    
end


surf(x,y,z,c);
colormap('jet');
if isfinite(customized_min(c,10))&isfinite(customized_max(c,10))&customized_max(c,10)>customized_min(c,10)
    caxis([customized_min(c,10) customized_max(c,10)]);
end
shading interp
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Altitude (m)');
hold on
tw=4;
xt=[0*ones(1,90);tw*ones(1,90);0*ones(1,90);-tw*ones(1,90)];
yt=[tw*ones(1,90);0*ones(1,90);-tw*ones(1,90);0*ones(1,90)];
zt=[1:90;1:90;1:90;1:90];
ct=customized_max(c,10)*[ones(1,90);ones(1,90);ones(1,90);ones(1,90)];
surf(xt,yt,zt,ct)

colorbar
% if isfinite(min(x))&isfinite(max(x))&isfinite(min(y))&isfinite(max(y))&isfinite(customized_min(z,1))&isfinite(customized_max(z,10))
    xlim([min(-tw,customized_min(x,0)) max(tw,customized_max(x,0))]);
    ylim([min(-tw,customized_min(y,0)) max(tw,customized_max(y,0))]);
if isfinite(customized_min(z,1))&isfinite(customized_max(z,10))    
    zlim([min(customized_min(z,1),minz) max(customized_max(z,10),maxz)]);
end

caxis([-2 2]);
view([-90,0]);