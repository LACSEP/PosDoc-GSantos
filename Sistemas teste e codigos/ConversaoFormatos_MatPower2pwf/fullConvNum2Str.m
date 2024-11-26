function [resConv] = fullConvNum2Str(t,x1,x2,c, fname)
nCarc1 = length(t);
nCarc2 = x2-x1+1;
resConv = [];
if nCarc2 > nCarc1
    difCar = nCarc2 - nCarc1;
    for k = 1 : 1 : difCar
        resConv = [resConv, ' '];
    end
    for k = 1 : 1 : nCarc1
        resConv = [resConv, t(k)];
    end
elseif nCarc2 == nCarc1
    resConv = t;
else
    fprintf('Excesso de caracteres na linha %d do cartao - campo: %s \n', c, fname);
    for k = 1 : 1 : nCarc2
        resConv = [resConv, t(k)];
    end
end