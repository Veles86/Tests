function [ phi, lambda ] = IG12toIGD12(E,N )

[phi_GRS80, lambda_GRS80,h_GRS80] = ITMtoGRS80(E,N);

[X_GRS80,Y_GRS80,Z_GRS80] = GG2XYZ(phi_GRS80,lambda_GRS80,h_GRS80,'GRS80');

[X_WGS84,Y_WGS84,Z_WGS84] = GRS80toWGS84(X_GRS80,Y_GRS80,Z_GRS80);

[phi_WGS84, lambda_WGS84] = XYZ2GG(X_WGS84,Y_WGS84,Z_WGS84,'WGS84');

phi = phi_WGS84;

lambda = lambda_WGS84;

end


function [ X,Y,Z ] = GG2XYZ( phi, lambda, h, datum )


phi = deg2rad(phi);
lambda = deg2rad(lambda);

if strcmp(datum,'WGS84')
    a = 6378137;
    f = 1 / 298.257223563;
elseif strcmp(datum,'GRS80')
    a = 6378137;
    f = 1/298.2572222;
else
    X = NaN;
    Y = NaN;
    Z = NaN;
    return;
end

sqr_e = 2 * f - f^2;
N = a ./ sqrt(1 - sqr_e * sin(phi).^2);

X = (N + h) .* cos(phi) .* cos(lambda);
Y = (N + h) .* cos(phi) .* sin(lambda);
Z = (N * (1 - sqr_e) + h) .* sin(phi);

end

function [X_WGS84,Y_WGS84,Z_WGS84] = GRS80toWGS84(X_GRS80,Y_GRS80,Z_GRS80)

dX=-24.002431;
dY=-17.103212;
dZ=-17.844441;
scale=1.000005424827;
RX=-0.000001600258;
RY=-0.000008982129;
RZ=0.000008094882;

% M = inv([scale RZ -RY; -RZ scale RX;RY -RX scale]);
% 
% X_WGS84 = M(1,1)*(X_GRS80-dX)+M(1,2)*(Y_GRS80-dY)+M(1,3)*(Z_GRS80-dZ);
% Y_WGS84 = M(2,1)*(X_GRS80-dX)+M(2,2)*(Y_GRS80-dY)+M(2,3)*(Z_GRS80-dZ);
% Z_WGS84 = M(3,1)*(X_GRS80-dX)+M(3,2)*(Y_GRS80-dY)+M(3,3)*(Z_GRS80-dZ);

M = ([scale RZ -RY; -RZ scale RX;RY -RX scale]);

Res = (M\[X_GRS80-dX,Y_GRS80-dY,Z_GRS80-dZ]')';
X_WGS84 = Res(:,1);
Y_WGS84 = Res(:,2);
Z_WGS84 = Res(:,3);

end

function [ X_GRS80,Y_GRS80,Z_GRS80 ] = WGS84toGRS80( X_WGS84,Y_WGS84,Z_WGS84 )

        dX=-24.002431;
        dY=-17.103212;
        dZ=-17.844441;
        scale=1.000005424827;
        RX=-0.000001600258;
        RY=-0.000008982129;
        RZ=0.000008094882;

        X_GRS80 = dX+(scale*X_WGS84+Y_WGS84*RZ-Z_WGS84*RY);
        Y_GRS80 = dY+(-X_WGS84*RZ+scale*Y_WGS84+Z_WGS84*RX);
        Z_GRS80 = dZ+(X_WGS84*RY-Y_WGS84*RX+scale*Z_WGS84);

end

function [ phi,lambda ] = XYZ2GG( X,Y,Z,datum )

if strcmp(datum,'WGS84')
    a = 6378137;
    f = 1 / 298.257223563;
elseif strcmp(datum,'GRS80')
    a = 6378137;
    f = 1/298.2572222;
else
    phi = NaN;
    lambda = NaN;
    return;
end


        lambda = atan(Y./X);

        sqr_e = 2*f-f^2;

        phi = atan((Z./sqrt(X.^2+Y.^2))*(1-sqr_e)^(-1));
        N = a ./ sqrt(1 - sqr_e * (sin(phi).^2));
        h = sqrt(X.^2 + Y.^2)./cos(phi)-N;

        stop = false;
        i=0;
        max_iterations = 100;

        while ~(stop)
        
            i=i+1;
            phi_prev = phi;
            h_prev = h;

            phi = atan((Z ./ sqrt(X.^2 + Y.^2)) .* (1 - sqr_e .* (N ./ (N + h))).^(-1));

            N = a ./ sqrt(1-sqr_e*(sin(phi).^2));
            h = sqrt(X.^2 + Y.^2)./cos(phi)-N;

            diff_phi = max(abs(phi - phi_prev));
            diff_h = max(abs(h - h_prev));

            if (((diff_phi<(1/206265))&&(diff_h<1e-4)) || i>max_iterations)
                stop = true;
            end
        end

phi = rad2deg(phi);
lambda = rad2deg(lambda);

end


function [ phi_GRS80,lambda_GRS80,h_GRS80 ] = ITMtoGRS80( y,x )

syms phi lambda

a = 6378137;
f = 1/298.2572222;

sqr_e = 2*f-f^2;

lambda0 = 0.614434732254689;


y0 = 219529.584;
x0 = 626907.390;

m0 = 1.0000067;

N = a./sqrt(1-(sqr_e)*sin(phi).^2);
J = (lambda - lambda0).*cos(phi);

t = tan(phi);
sqr_ni = (sqr_e*cos(phi).^2)/(1-sqr_e);

D2 = (J.^2).*(1-t.^2+sqr_ni)/6;
D3 = (J.^4).*(5-18*t.^2+14*sqr_ni+t.^4-58*(t.^2).*sqr_ni)/120;

C2 = (J.^2).*(5-t.^2+9*sqr_ni+4*sqr_ni.^2)/12;
C3 = (J.^4).*(61-58*t.^2+270*sqr_ni+t.^4-330*(t.^2).*sqr_ni)/360;

CA = 1.005052500;
CB = 0.002531553;
CC = 0.000002657;
CD = 0.000000003;


Sm0 = 3512424.3388;


Sm = m0*a*(1-sqr_e)*(CA*phi-CB*sin(2*phi)+CC*sin(4 * phi)-CD*sin(6 * phi));



temp_phi = zeros(length(y),1);
temp_lambda = zeros(length(y),1);
for i = 1: length(y)
    eqns = [y(i) == y0 + m0*N.*J.*(1+D2+D3);x(i) == x0 + (Sm-Sm0)+m0.*N.*(J.^2).*t.*(1+C2+C3)/2];
S = vpasolve(eqns);
temp_phi(i) = double(S.phi);
temp_lambda(i) = double(S.lambda);

end
phi_GRS80 = rad2deg(temp_phi);
lambda_GRS80 = rad2deg(temp_lambda);
h_GRS80=zeros(size(phi_GRS80,1),1);


end




function [ y, x ] = GRS80toITM( phi,lambda )

        a = 6378137;
        f = 1/298.2572222;
        
        sqr_e = 2*f-f^2;

        lambda0 = 0.614434732254689;


        y0 = 219529.584;
        x0 = 626907.390;

        m0 = 1.0000067;

        N = a./sqrt(1-(sqr_e)*sin(phi).^2);
        J = (lambda - lambda0).*cos(phi);

        t = tan(phi);
        sqr_ni = (sqr_e*cos(phi).^2)/(1-sqr_e);

        D2 = (J.^2).*(1-t.^2+sqr_ni)/6;
        D3 = (J.^4).*(5-18*t.^2+14*sqr_ni+t.^4-58*(t.^2).*sqr_ni)/120;

        C2 = (J.^2).*(5-t.^2+9*sqr_ni+4*sqr_ni.^2)/12;
        C3 = (J.^4).*(61-58*t.^2+270*sqr_ni+t.^4-330*(t.^2).*sqr_ni)/360;

        CA = 1.005052500;
        CB = 0.002531553;
        CC = 0.000002657;
        CD = 0.000000003;


        Sm0 = 3512424.3388;


        Sm = m0*a*(1-sqr_e)*(CA*phi-CB*sin(2*phi)+CC*sin(4 * phi)-CD*sin(6 * phi));

        y = y0 + m0*N.*J.*(1+D2+D3);
        x = x0 + (Sm-Sm0)+m0.*N.*(J.^2).*t.*(1+C2+C3)/2;

end





