function xf = forcaBruta
    x = 17*ones(1,345);  
    xl = x;
    for i = 1:length(x)
        right = 0;
        xl(i) = x(i) - 1;
        while right == 0
            [~,csig,cu] = ISCSO_2021(xl,0);
            if csig == 0 && cu == 0
                if xl(i) > 1
                    x(i)  = x(i) - 1;
                    xl(i) = x(i) - 1;
                else
                    right = 1;
                    xf(i) = xl(i);
                end
            else
                xf(i) = x(i);
                xl(i) = x(i);
                right = 1;
            end  
        end
    end

end