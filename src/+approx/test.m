classdef test < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function randomSeed(testCase)
            rng(1);
        end
    end
    
    methods (Test)
        function chebvander(testCase)
            [~, V] = approx.chebvander(10);
            testCase.verifyEqual(cond(V), sqrt(2), 'abstol', eps(100));
        end
        
        function chebInterp1D(testCase)
            p = 30;
            [rvec, V] = approx.chebvander(p);
            x = rand(100, 1)*2 - 1;
            E = approx.chebevalmat(x, p);
            f = @(x) exp(-(x+1).^2);
            finterp = E*(V\f(rvec));
            testCase.verifyEqual(finterp, f(x), 'reltol', eps(100));
        end
        
        function kronMatApply2Dsame(testCase)
            p = 10;
            A = rand(p, p);
            x = rand(p^2, 1);
            y1 = kron(A, A)*x;
            y2 = approx.kronmat_apply(A, x, 2);
            testCase.verifyEqual(y1, y2, 'reltol', eps(p));
        end

        function kronMatApply2D(testCase)
            m1 = 10;
            m2 = 11;
            n1 = 12;
            n2 = 13;
            A = rand(m1, n1);
            B = rand(m2, n2);
            x = rand(n1*n2, 1);
            y1 = kron(A, B)*x;
            y2 = approx.kronmat2_apply(A, B, x);
            testCase.verifyEqual(y1, y2, 'reltol', eps(10));
        end

        function kronMatApply3Dsame(testCase)
            p = 10;
            A = rand(p, p);
            x = rand(p^3, 1);
            y1 = kron(A, kron(A, A))*x;
            y2 = approx.kronmat_apply(A, x, 3);
            testCase.verifyEqual(y1, y2, 'reltol', eps(p));
        end

        function kronMatApply3D(testCase)
            m1 = 3;
            m2 = 5;
            m3 = 7;
            n1 = 8;
            n2 = 13;
            n3 = 2;
            A = rand(m1, n1);
            B = rand(m2, n2);
            C = rand(m3, n3);
            x = rand(n1*n2*n3, 1);
            M = kron(kron(A, B), C);
            y1 = M*x;
            y2 = approx.kronmat3_apply(A, B, C, x);
            testCase.verifyEqual(y1, y2, 'reltol', eps(10));
        end
        
        function kronMatApply4Dsame(testCase)
            p = 5;
            A = rand(p, p);
            x = rand(p^4, 1);
            y1 = kron(A, kron(A, kron(A, A)))*x;
            y2 = approx.kronmat_apply(A, x, 4);
            testCase.verifyEqual(y1, y2, 'reltol', eps(p));
        end

        function chebevalmat3_fast(testCase)
            p = 30;
            N = 100;
            M = rand(p, p);
            x = rand(N, 1)*2 - 1;
            y = rand(N, 1)*2 - 1;
            z = rand(N, 1)*2 - 1;
            E = approx.chebevalmat3(x, y, z, p, M);
            ET = transpose(E);
            f = rand(100, 1);
            g1 = ET*f;
            % Trans apply with weight
            g2 = approx.chebevalmat3_trans_apply(x, y, z, p, f, M);
            testCase.verifyEqual(g1, g2, 'abstol', eps(p*N));
            expa = rand(p^3, 1);
            % Straight apply without weight
            E1 = approx.chebevalmat3(x, y, z, p);
            fi = approx.chebevalmat3_apply(x, y, z, p, expa);
            testCase.verifyEqual(fi, E1*expa, 'abstol', eps(p*N));
            % Test that dimensions get right with scalar eval
            fi1 = approx.chebevalmat3_apply(x(1), y(1), z(1), p, expa);
            testCase.verifyEqual(fi1, fi(1), 'abstol', eps(10));
        end
        
        function chebInterp2D(testCase)
            p = 30;
            N = 100;
            x = rand(N, 1)*2 - 1;
            y = rand(N, 1)*2 - 1;
            f = @(x, y) exp(-(x + y/2).^2);
            [rvec, V] = approx.chebvander(p);
            [xp, yp] = ndgrid(rvec, rvec);
            fp = f(xp, yp);
            Vi = inv(V);
            E = approx.chebevalmat2(x, y, p);
            fi = E * approx.kronmat_apply(Vi, fp(:), 2);
            testCase.verifyEqual(fi, f(x, y), 'reltol', eps(1000));
        end

        function chebInterp3D(testCase)
        % Test tenstor product cube -> NU points
            p = 30;
            N = 100;
            x = rand(N, 1)*2 - 1;
            y = rand(N, 1)*2 - 1;
            z = rand(N, 1)*2 - 1;
            f = @(x, y, z) exp(-(x + y/2 + z/3).^2);
            [rvec, V] = approx.chebvander(p);
            [xp, yp, zp] = ndgrid(rvec, rvec, rvec);
            fp = f(xp, yp, zp);
            Vi = inv(V);
            E = approx.chebevalmat3(x, y, z, p);
            expa = approx.kronmat_apply(Vi, fp(:), 3);
            fi = E * expa;
            testCase.verifyEqual(fi, f(x, y, z), 'reltol', eps(1000));
            % Fast apply
            fi2 = approx.chebevalmat3_apply(x, y, z, p, expa);
            testCase.verifyEqual(fi, fi2, 'abstol', eps(10));
        end
        
        function chebInterp3Dgrid2grid(testCase)
        % Test tensor product cube -> tensor product rectangle
            p = 25;
            q = 5;
            f = @(x, y, z) exp(-(x + y/2 + z/3).^2);
            [rvec, V] = approx.chebvander(p);
            [xp, yp, zp] = ndgrid(rvec, rvec, rvec);
            fp = f(xp, yp, zp);
            xq = (chebpts(q, 1)+1)/2;
            yq = (chebpts(q+1, 1)-1)/2;
            zq = chebpts(q+2, 1)/2;
            Ex = approx.chebevalmat(xq, p);
            Ey = approx.chebevalmat(yq, p);
            Ez = approx.chebevalmat(zq, p);
            [x, y, z] = ndgrid(xq, yq, zq);
            Ux = Ex/V;
            Uy = Ey/V;
            Uz = Ez/V;
            fi = approx.kronmat3_apply(Uz, Uy, Ux, fp);
            fref = f(x, y, z);
            testCase.verifyEqual(fi, fref(:), 'reltol', eps(1000));
        end


        
    end
end
   
