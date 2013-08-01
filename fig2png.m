clear all;close all;
ad='/home/schmmark/superposed/';%input address to sequences
save = '/home/schmmark/superposed/images/';

files=dir(fullfile(ad,'*.fig') );
num_file=size(files);
frames_cnt=1;

for i=1:num_file
    open([ad files(i).name]);
    set(gcf,'PaperPositionMode','auto')
    print('-dpng',...   % print png file
    fullfile(save,sprintf('%s%s',files(i).name,'.png')))%,'-zbuffer','-r200');
    close all
end