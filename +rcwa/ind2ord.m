function [ord] = ind2ord(ind, sys)
    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    if ~isempty(find(ind > sys.nTot))
        error('RCWA:ind2ord', 'Requested Orders excede those calculated')
    end

    ind = int32(ind);

    NM = 2 * int32(sys.maxM) + 1;

    x = mod(ind - 1, NM(1)) - sys.maxM(1);
    y = idivide(ind - 1, NM(1)) - sys.maxM(2);
    ord = [x; y]';
end
