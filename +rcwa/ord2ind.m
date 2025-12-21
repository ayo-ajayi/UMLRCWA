function ind = ord2ind(ord, sys)

    % converts diffraction order pairs (m,n) to linear index, i.e for mapping function between two ways of labelling diffraction orders
    % (m,n)=diffraction orders in x and y directions
    % for example: maxM=[2,2]:
    % (-2,-2)  (-1,-2)  (0,-2)  (1,-2)  (2,-2)
    % (-2,-1)  (-1,-1)  (0,-1)  (1,-1)  (2,-1)
    % (-2, 0)  (-1, 0)  (0, 0)  (1, 0)  (2, 0)  ← (0,0) is specular
    % (-2, 1)  (-1, 1)  (0, 1)  (1, 1)  (2, 1)
    % (-2, 2)  (-1, 2)  (0, 2)  (1, 2)  (2, 2)
    % in linear index becomes:
    %         1    2    3    4    5
    %         6    7    8    9   10
    %        11   12   13   14   15  ← Index 13 is (0,0)
    %        16   17   18   19   20
    %        21   22   23   24   25

    % Copyright (c) 2014-2017 by University of Massachusetts Lowell

    if sum(ord(:, 1) > sys.maxM(1) | ord(:, 2) > sys.maxM(2))
        error('RCWA:ind2ord', 'Requested Orders excede those calculated')
    end

    aa = ord(:, 1) + sys.maxM(1) + 1;
    bb = ord(:, 2) + sys.maxM(2);

    ind = aa + bb * (2 * sys.maxM(1) + 1);
    % returns a linear index in the form (N x 1) vector
end
