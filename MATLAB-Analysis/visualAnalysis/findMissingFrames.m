a = ones(size(m.pdRise));
b = 2*ones(size(m.pdFall));
c = [a b];
[~, temp] = sort([m.pdRise m.pdFall]);
c = c(temp);
discontLoc = find(c == 0);