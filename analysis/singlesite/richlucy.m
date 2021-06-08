function out = richlucy(in,psf,it)
%   out = richlucy(in,psf,it)
%
% Iterative Richardson Lucy sharpening algorithm. Sharpens an image 'in'
% through 'it' iterations of the Richardson-Lucy algorithm with the
% psf/kernel 'psf'. The algorithm is in theoretically norm preserving;
% boundary effects may, however, violate this in the present
% implementation. [ST 2015-10]
%
% Compare to Matlab's deconvlucy.m (imaging toolbox)

    out = in;

    for j = 1:it
        out = out.*(conv2(in./conv2(out,psf,'same'),psf,'same'));
    end


end