function [n, k] = eps2nk(epsRe, epsIm)
    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    if nargin == 1
        er = imag(epsRe);
        ei = real(epsIm);
    elseif nargin == 2
        er = epsRe;
        ei = epsIm;
    end

    x = sqrt(er.^2 + ei.^2);
    n = sqrt((x + er) / 2);
    k = sqrt((x - er) / 2);

    if nargout == 1
        n = n + 1i * k;
    end

end
