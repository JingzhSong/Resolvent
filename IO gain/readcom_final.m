function [Cm,Cd,Ck,Y] = readcom_final(method,eddy,wall,selc,effect,kx,omega)

% method = 'IOA' or 'SIOA';
% eddy = 'eddyon' or 'eddyoff';
% wall = 'com';
% selc = 'RS' or 'SIN';
% effect = 'good' or 'bad';

if strcmp(wall,'com')
    if strcmp(effect,'good')
        if strcmp(selc, 'SIN')
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=-5.050505e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=-3.030303e-01;Ck=2.885859e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 12
                % Cd=-3.737374e+00;Ck=2.885859e+04;
                Cd=-4.141414e+00;Ck=2.893939e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 12
                % Cd=7.070707e-01;Ck=2.889899e+04;
                Cd=9.090909e-01;Ck=2.893939e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=5.050505e-01;Ck=5.090909e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=7.070707e-01;Ck=5.090909e+02;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 1
                Cd=-8.181818e+00;Ck=4.101010e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 1
                % Cd=-3.030303e-01;Ck=4.949495e+02;
                Cd=5.555556e+00;Ck=4.914141e+02;
            end
        end
        
        if strcmp(selc, 'RS')
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=-3.030303e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=-9.090909e-01;Ck=2.845455e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 12
                Cd=-1.313131e+00;Ck=2.889899e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 12
                  Cd=-1.515152e+00;Ck=2.845455e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=3.030303e-01;Ck=5.126263e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=-1.111111e+00;Ck=4.136364e+02;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 1
                Cd=-1.010101e-01;Ck=5.090909e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 1
                Cd=-2.121212e+00;Ck=4.207071e+02;
            end
        end

    end
    
    if strcmp(effect,'bad')
                if strcmp(selc, 'SIN')
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=1.010101e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=3.030303e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 12
                Cd=-1.010101e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 12
                Cd=-1.010101e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=-5.050505e-01;Ck=5.161616e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 1
               Cd=-5.050505e-01;Ck=5.161616e+02;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 1
                Cd=3.030303e-01;Ck=5.161616e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 1
                Cd=-3.131313e+00;Ck=5.161616e+02;
            end
        end
        
        if strcmp(selc, 'RS')
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 12
                Cd=3.030303e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 12
%                Cd=-9.090909e-01;Ck=2.845455e+04;
               Cd=-7.070707e-01;Ck=2.841414e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 12
                Cd=-1.010101e-01;Ck=2.881818e+04;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 12
                Cd=-1.313131e+00;Ck=2.845455e+04;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=-5.050505e-01;Ck=5.161616e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyoff') && kx == 1
                Cd=9.090909e-01;Ck=4.101010e+02;
            end
            if strcmp(method, 'IOA') && strcmp(eddy, 'eddyon') && kx == 1
               Cd=3.030303e-01;Ck=5.161616e+02;
            end
            if strcmp(method, 'SIOA') && strcmp(eddy, 'eddyon') && kx == 1
%                 Cd=-2.121212e+00;Ck=4.207071e+02;
               Cd=-1.010101e-01;Ck=4.136364e+02;
            end
        end
    end
    Cm = 2;
    
    [ReY,ImY] = C2Y(Cm,Ck,Cd,omega);
    fprintf('ReY=%d\nImY=%d',ReY,ImY);
    
    Y = ReY + 1i*ImY;
else
    Y = [];
    Cm = [];
    Cd = [];
    Ck = [];
end
end

