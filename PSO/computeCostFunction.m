function [f,c1,c2] = computeCostFunction(x)
    x = x*(37-1) + 1;
    x = round(x);
    x = max(1,min(37,x));
    [f,c1,c2] = ISCSO_2021(x,0);
%     v1 = 0;
%     v2 = 0;
%     if c1 > 0
%         v1 = 1;
%     end
%     if c2 > 0
%         v2 = 1;
%     end
    f = f + (c1 + c2)*1e10;
end