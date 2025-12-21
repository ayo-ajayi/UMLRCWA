function [epsRe, epsIm] = nk2eps(nI, kI)
    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    if nargin == 1
        k = imag(nI);
        n = real(nI);
    elseif nargin == 2
        n = nI;
        k = kI;
    end

    epsRe = n.^2 - k.^2;
    epsIm = 2 * n .* k;

    if nargout == 1 || nargout == 0
        epsRe = epsRe + 1i * epsIm;
    end

end
