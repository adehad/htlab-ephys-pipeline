m.pdRise(2,:) = 1;
m.pdFall(2,:) = 2;
interlaced = zeros(2,max([length(m.pdRise) length(m.pdFall)])*2);
interlaced(:,1:2:end-2) = m.pdRise;
interlaced(:,2:2:end) = m.pdFall;
temp = diff(interlaced(2,:));
discontLoc = find(temp == 0)