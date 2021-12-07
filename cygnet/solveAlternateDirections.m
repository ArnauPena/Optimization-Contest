%% Alternate  Directions method
% alternes. primer proves una s,busques un tau tal que canviant la s amb la
% derivada i fent que sigui entre 0 i 1 la funcio objectiu baixi
% un cop trobat, busques amb la nova s, mirem restriccions. es compleixen?
% si, perfecte. no? busques una s tal que projectes i puges una mica la s
% amb el gradient pq compleixi la restriccio, de manera diferent a
% l'augmented

function solveAlternateDirections
    [sec]   = computeSectionData;
    tol     = 1e-3;
    s0      = 0.99*ones(1,345);
    [c0,c_u0,c_sig0] = computeParameters(s0);

    tau1 = 0.5;
    tau2 = 0.5;

    rightTau1 = 0;
    dSTau = 1;
    while dSTau > tol
        while rightTau1 == 0
    %         smean = s0 - tau1 * dC0;
            [dA0, ~, ~] = computeParametersGradient(s0,sec);
            smean = s0 - tau1 * dA0;
            sTau1 = max (0, max(0, min(1,smean)));
            x = calculateSectionID(sTau1);
            [c,~,~] = ISCSO_2021(x,0);
            if c < c0
                rightTau1 = 1;
            else
                tau1 = tau1 + 0.1;
            end
        end
    
        rightTau2 = 0;
        while rightTau2 == 0
            A0 = computeSection(s0,sec);
            [dA0, dC_u0, dC_sig0] = computeParametersGradient(s0,sec);
            smean = s0 - tau2 .* (dA0 ./ (A0).^2);
            sTau2 = max (0, max(0, min(1,smean)));
            [~,c_u,c_sig] = computeParameters(sTau2); % entenc que es stau2
    
    %         if c_u < c_u0 && c_sig < c_sig0
            % Pot ser que s'hi hagi d'afegir l'igual?
            if c_u <= 0 && c_sig <= 0
                rightTau2 = 1;
            else
                tau2 = tau2 + 0.1;
            end
        end
    
        dSTau = norm(sTau1-sTau2)/norm(sTau1);
        rightTau1 = 0;
        rightTau2 = 0;
        s0 = sTau1;
    end

end

function x = calculateSectionID(s)
    x = s*(37-1) + 1;
    x = round(x);
end

function [c,c_u,c_sig] = computeParameters(s)
    x = calculateSectionID(s);
    [w,v1,v2] = ISCSO_2021(x,0);
    c     = w;
    c_sig = v1;
    c_u   = v2;       
end

function [dc,dc_u,dc_sig] = computeParametersGradient(s,sec)
    dA     = computeSectionGradient(s,sec);
    A      = computeSection(s,sec);
    dc     = dA;
    dc_u   = -dA./A.^2;
    dc_sig = -dA./A.^2;
end

function [Si] = computeData()
    run("DataSections.m")
end

function [sections] = computeSectionData()
    Si = computeData();
    sections = Si(:,1);
end

function A = computeSection(s,sec)
    sections = sec;
    Amin     = min(sections);
    Amax     = max(sections);
    p        = 3;
    A        = Amin * (1-s.^p) + s.^p * Amax;
end

function gradA = computeSectionGradient(s,sec)
    sections = sec;
    Amin     = min(sections);
    Amax     = max(sections);
    p        = 3;
    gradA    = - Amin * p * s.^(p-1) + Amax * p*s.^(p-1);
end

function [dC, dV1, dV2] = calculateFiniteGradient(s)
    x0 = calculateSectionID(s);
    [c0, v10, v20] = ISCSO_2021(x0,0);
    for i = 1:numel(s)
        s(i) = s(i)+1/37;
        x = calculateSectionID(s);
        [c, v1, v2] = ISCSO_2021(x,0);
        ds = 1/37;
        dC(i)  = (c-c0) / (ds);
        dV1(i) = (v1-v10) / (ds);
        dV2(i) = (v2-v20) / (ds);
    end
end

function [cost, const_u, const_sig, lu, lsig, tauV] = computeMonitoring(c, c_u, c_sig, l_u, l_sig, tau, w0, iter, s,...
    cost, const_u, const_sig, lu, lsig, tauV, monitor)
        cost(iter+1)      = c;
        const_u(iter+1)   = c_u;
        const_sig(iter+1) = c_sig;
        lu(iter+1)        = l_u;
        lsig(iter+1)      = l_sig;
        tauV(iter+1)      = tau;
        x                 = calculateSectionID(s);
        auglagr_monitoring(monitor, 1:(iter+1), x, cost, const_sig, const_u, lu, lsig ,tauV, w0);
end
