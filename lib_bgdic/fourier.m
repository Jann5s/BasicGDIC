function [f, yfit] = fourier(x,y,n)
% [f, yfit] = fourier(x,y,n)
%
% fourier transform of the signal (x,y) using n wavelengths
%
% f is a strucutre contiaining
% f.y0, the dc component
% f.a0, the zero order wave
% f.a, f.b, the amplitudes (cos, sin, respectively) of each wavelength
%
% yfit returns the fourier fit of the data, i.e.
%     for k = 1:n
%         yfit = (a0 / 2) + y0;
%         for k = 1:n
%             yfit = yfit + a(k)*cos(k*x) + b(k)*sin(k*x);
%         end
%     end

% start at zero
y0 = y(1);
y  = y - y0 ;

% integration spacing
dx = mean(diff(x));

% get the fourier parameters
a0 = (1/pi) * sum(y) * dx;
a = zeros(1,n);
b = zeros(1,n);
for k = 1:n
    a(k) = (1/pi) * sum( y .* cos(k*x) ) * dx;
    b(k) = (1/pi) * sum( y .* sin(k*x) ) * dx;
    
end

f.y0 = y0;
f.a0 = a0;
f.a = a;
f.b = b;

if nargout == 1
    return
elseif nargout == 2
    % Reconstruction
    for k = 1:n
        yfit = (a0 / 2) + y0;
        for k = 1:n
            yfit = yfit + a(k)*cos(k*x) + b(k)*sin(k*x);
        end
    end
end
