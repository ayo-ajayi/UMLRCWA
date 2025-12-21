function out = mySqrt(eps, branchcut)
    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    % calculates sqrt with a specified branchcut
    if nargin == 2
        temp = mod(angle(eps) - branchcut, 2 * pi) + branchcut;
        out = sqrt(abs(eps)) .* exp(temp / 2 * 1i);
    else
        out = sqrt(eps);
    end

end
