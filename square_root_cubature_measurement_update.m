function [xS_corrected, M_corrected] = square_root_cubature_measurement_update(xS0,measurement,measurement_fun,measurement_noise)

%Performs the measurement update of the CD-CKF in square-root form based on
%impmentation of I. Arasaratnam.
%   Input Arguments:
%   xS0: average and a square root of predicted covariance matrix of class
%   mean_covariance_sqrt_cls
%   measurement: measurement at step
%   measurement_fun: (non-stochastic) measurement function z = f(x)
%   measurement_noise: Covariance matrix of measurement noise (not a square
%   root)
%   Output Argument:
%   xS_corrected: corrected average and square root covariance matrix of
%   class mean_covariance_sqrt_cls

%Translate variable names and formats:

Rsqrt = chol(measurement_noise,'lower');
% nx = xS0.dim;
% xkk1 = xS0.mean;
% Skk1 = xS0.c_sqrt;

nx = size(xS0,1);
xkk1 = xS0(:,1);
Skk1 = xS0(:,2:4);

nz = numel(measurement_fun(xkk1));

nPts = 2 * nx;
QPtArray = sqrt(nx) * cat(2,eye(nx),-eye(nx));

Xi =  repmat(xkk1,1,nPts) + Skk1*QPtArray;
%Finds quadrature points.
Zi = zeros(nz,nPts);
for j = 1:nPts
    Zi(:,j) = measurement_fun(Xi(:,j));
end
%Finds expected measurement (Update: accomodates measurement function that
%only takes a single state variable.)
zkk1 = sum(Zi,2)/nPts;  %%% zkk1 represents y_bar in eq.48 in manuscript
%Finds averaged measurement
X = (Xi-xkk1)/sqrt(nPts);
%X = sqrt(2) * [Skk1| -Skk1]
Z = (Zi-zkk1)/sqrt(nPts);  
%Z gives corresponding measurements
[~,S] = qr([Z Rsqrt; X zeros(nx,nz)]',0);
%Uses QR factorization to drop singular dimensions
S = S';
%Takes "lower triangular" form due to convention of square root
X = S(1:nz,1:nz);

Y = S(nz+1:end,1:nz);

Z = S(nz+1:end,nz+1:end);
%Note: not of same dimension as the previous Z matrix.
G = Y/X;

% Skk = Z;
xkk = xkk1 + G*(measurement-zkk1);  
                      
%Reformat output:
xS_corrected = xkk;
M_corrected = Z;
           
end
