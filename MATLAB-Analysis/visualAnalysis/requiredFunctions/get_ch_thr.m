function thr = get_ch_thr(fname,nch,tch,ndata,chname,f,thr,bgoff,recenter)
%
%allow user to set event detection threshold for a channel
%
%IN: fname = file to load
%    nch = total ch in file
%    tch = channel to set threshold
%    ndata = data load block size (pts)
%    chname = string name of ch
%    f = denoising filter for ch
%    thr = initial thr for ch
%    bgoff = subtract bg (1=yes, 0=no)
%
%OUT: thr = final thr for ch
%
%Requires: get_events.m

flen = length(f);
fid = fopen(fname,'r');
cmdlist = sprintf('channel: %s\n\n(s)et new threshold\n(l)oad next data block\n(a)ccept settings\nset (b)lock size',chname);
cmd = 'l';
while cmd ~= 'a'
    switch cmd
        case 'l'
            x = fread(fid,[nch ndata],'int16');
            if bgoff
                xmu = mean(x(tch,:));
            else xmu = 0;
            end
            [events,waves,xi] = get_events(x,tch,f,10,10,thr,bgoff,recenter);
            clf;
            hold on
            plot(x(tch,:)-xmu,'b')
            plot(xi,'r')
            plot(events,xi(events),'m*');
            plot([0 ndata],[1 1]*thr,'k:')
            zoom on;
            title(sprintf('%s\nfound %d events',cmdlist,length(events)));
        case 's'
            thr = input(sprintf('Enter new threshold (old= %d): ',thr));
            [events,waves,xi] = get_events(x,tch,f,10,10,thr,bgoff,recenter);
            clf;
            hold on
            plot(x(tch,:)-xmu,'b')
            plot(xi,'r')
            plot(events,xi(events),'m*');
            plot([0 ndata],[1 1]*thr,'k:')
            zoom on;
            title(sprintf('%s\nfound %d events',cmdlist,length(events)));
        case 'b'
            ndata = input(sprintf('Enter new data load block size (old= %d): ',ndata));
    end
    %     [x1,y1,cmd] = ginput2(1);
    cmd = input('next command: ','s');
end
fclose(fid);
