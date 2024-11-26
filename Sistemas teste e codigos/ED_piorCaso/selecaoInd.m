function [sel] = selecaoInd(nInd, targetInd)

aux = [1:nInd];
aux(targetInd) = [];

for i = 1 : 1 : 5
    k = randi ([1 nInd-i]);
    pos = aux(k);
    aux(k) = [];
    sel(i) = pos;
end

