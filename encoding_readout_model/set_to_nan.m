function [a,b]=set_to_nan(a,b)

if isnan(a) || isnan(b)
    a=nan;
    b=nan;
end

end