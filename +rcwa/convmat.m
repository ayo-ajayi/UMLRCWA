% Copyright (c) 2014-2017 by University of Massachusetts Lowell
function [C, Cinv] = convmat(A, maxM, debug)
    %Given a 1D/2D spacial representation of the data
    %Perform 1D/2D fourier transform
    %Returns convolution matrix
    %Optional returns the inverse convolution matrix.
    %Currently only works for 2D.

    NM = 2 * maxM + 1;
    nModeTot = prod(NM);

    if length(NM) == 1
        nx = NM;
        ny = 1;
    elseif length(NM) == 2
        nx = NM(1);
        ny = NM(2);
    end

    x = 1:nx * ny;

    indX = mod(x - 1, nx) + nx;
    indY = ceil(x / nx) + ny;

    sA = size(A);
    i0 = 1 + floor(sA / 2);

    %The following command is strange but dut to conventions in the
    %cooley-turkey algorithm, and in particular the fftw library + matlab, this
    %is the corrent way to take an fft of an image centered at (0,0), AND match
    %the expectations for the analytic k-space fouier transform we are used to.

    %Calculate the fourier convolution matrix of A(x,y)
    if length(NM) == 2
        FA = ifftshift(ifft2(ifftshift(A)));
    elseif length(NM) == 1
        FA = ifftshift(ifft(ifftshift(A)));
    end

    C = FA(bsxfun(@minus, indY', indY) + i0(1) + (bsxfun(@minus, indX', indX) + i0(2) - 1) * sA(1));

    %Calculate the fourier convolution matrix of 1/A(x,y)
    if nargout == 2

        if length(NM) == 2
            FA = ifftshift(ifft2(ifftshift(1 ./ A)));
        elseif length(NM) == 1
            FA = ifftshift(ifft(ifftshift(1 ./ A)));
        end

        Cinv = FA(bsxfun(@minus, indY', indY) + i0(1) + (bsxfun(@minus, indX', indX) + i0(2) - 1) * sA(1));
    end

    %%%%%%%%%%%% Debuging Code %%%%%%%%%%%%%%%%%%
    if (nargin == 3 && debug == 1)
        xA=-sA(1) / 2:sA(1) / 2 - 1;
        yA=-sA(2) / 2:sA(2) / 2 - 1;
        figure()
        subplot(2, 3, 1)
        imagesc(xA, yA, real(FA))
        xlim([-2 * maxM(1), 2 * maxM(1)])
        ylim([-2 * maxM(2), 2 * maxM(2)])
        shading flat
        colorbar
        subplot(2, 3, 2)
        imagesc(xA, yA, imag(FA))
        xlim([-2 * maxM(1), 2 * maxM(1)])
        ylim([-2 * maxM(2), 2 * maxM(2)])
        shading flat
        colorbar

        subplot(2, 3, 3)
        imagesc(real(C))
        shading flat
        colorbar
        subplot(2, 3, 4)
        imagesc(imag(C))
        shading flat
        colorbar

        if nargout == 2
            subplot(2, 3, 5)
            imagesc(real(Cinv))
            shading flat
            colorbar
            subplot(2, 3, 6)
            imagesc(imag(Cinv))
            shading flat
            colorbar
        end

    end

end
