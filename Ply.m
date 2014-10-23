classdef Ply
    properties
        Material
        Theta
    end
    methods
        % Main constructor
        function p = Ply(nMaterial, nTheta)
            p.Material = nMaterial;
            p.Theta = nTheta;
        end    
        
        % Check if input stress vector leads to ply failure according to
        % one of FailureCriterion
        % StressVector is [sigmaX, sigmaY, tauXY]
        % Possible criterions: MaximumStrength, MaximumStrain, Tsai-Hill, Tsai-Wu
        function out = CheckFailure(obj, StressVector, FailureCriterion)
            sigmaX1 = StressVector(1);
            sigmaY1 = StressVector(2);
            tauXY1 = StressVector(3);
            
            switch FailureCriterion
                case 'MaximumStrength'
                    isFailure = @(sigmaX, sigmaY, tauXY)maxStrength(sigmaX, sigmaY, tauXY, obj);

                case 'MaximumStrain'
                    isFailure = @(sigmaX, sigmaY, tauXY)maxStrain(sigmaX, sigmaY, tauXY, obj);

                case 'Tsai-Hill'
                    isFailure = @(sigmaX, sigmaY, tauXY)TsaiHill(sigmaX, sigmaY, tauXY, obj);


                case 'Tsai-Wu'
                    isFailure = @(sigmaX, sigmaY, tauXY)TsaiWu(sigmaX, sigmaY, tauXY, obj);

                otherwise
                    error([FailureCriterion ' is not a valid failure criterion. Possible criterions: MaximumStrength, MaximumStrain, Tsai-Hill, Tsai-Wu'])
            end 
            
            out = isFailure(sigmaX1, sigmaY1, tauXY1);
        end
        
        % Draw failure envelope for given tauXY
        % and failure criterion
        % Possible criterions: MaximumStrength, MaximumStrain, Tsai-Hill, Tsai-Wu
        function [envelope, x, y, tau] = DrawFailureEnvelope(obj, tau, FailureCriterion)
            x = -1e3:20e0:1e3;
            y = x;
            envelope = zeros (length(x), length(y));
            
            NTotal = length(x)*length(y);
            k = 0;
            progressbar('Building failure envelope...');
            for i=1:length(x)
                for j = 1:length(y)
                    envelope(i, j) = obj.CheckFailure([x(i), y(j), tau], FailureCriterion);
                    k = k + 1;
                    progressbar(k/NTotal);
                end
            end
            progressbar(1);

            envelope(envelope == true) = NaN;

            f = figure;
            pc = pcolor(x, y, envelope);
            set(pc,'LineStyle','none');
            set(gca, 'XMinorTick','on','FontWeight','bold',...
                'FontSize',14,'YGrid','on','XGrid','on');
            ylabel('\sigma_y');
            xlabel('\sigma_x');
            title({...
                    ['Failure envelope for ' FailureCriterion ' criterion'];... 
                    ['Single-ply ' obj.Material.Name ', \theta = ' num2str(obj.Theta*180.0/pi, '%.1f') char(176)]});
            %saveas(f, criterion, 'png');
        end
        
        % Calculates transformation matrix for off-axis ply
        function T = TransformationMatrix(obj)
            c = cos(obj.Theta);
            s = sin(obj.Theta);

            T = [...
                    c*c,    s*s,    2*s*c;...
                    s*s,    c*c,    -2*s*c;...
                    -s*c,   s*c,    c*c - s*s];
        end
        
        % Calculates transformed reduced stiffness matrix
        function out = QBar(obj)
            theta = obj.Theta;
            
            v21 = obj.Material.E2 * obj.Material.v12 / obj.Material.E1;
            Q11 = obj.Material.E1 / (1 - obj.Material.v12 * v21);
            Q22 = obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q12 = obj.Material.v12 * obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q66 = obj.Material.G12;

            U1 = (3*Q11 + 3*Q22 + 2*Q12 + 4*Q66) / 8;
            U2 = (Q11 - Q22) / 2;
            U3 = (Q11 + Q22 - 2*Q12 - 4*Q66) / 8;
            U4 = (Q11 + Q22 + 6*Q12 - 4*Q66) / 8;
            U5 = (Q11 + Q22 - 2*Q12 + 4*Q66) / 8;

            Q11dash = U1 + U2*cos(2*theta) + U3*cos(4*theta);
            Q12dash = U4 - U3*cos(4*theta);
            Q22dash = U1 - U2*cos(2*theta) + U3*cos(4*theta);
            Q66dash = U5 - U3*cos(4*theta);
            Q16dash = 0.5*U2*sin(2*theta) + U3*sin(4*theta);
            Q26dash = 0.5*U2*sin(2*theta) - U3*sin(4*theta);

            out = [Q11dash Q12dash Q16dash; Q12dash Q22dash Q26dash; Q16dash Q26dash Q66dash];
        end
        
        % Calculates strains from provided stresses and ply properties
        function strains = GetStrainsFromStress(obj, StressVector)
            qb = obj.QBar();
            strains = qb \ reshape(StressVector, [3,1]); % The same as inv(qb) * StressVector
        end
    end
end







% Helper methods

function isFailure = maxStrength(sigmaX, sigmaY, tauXY, ply)
    T = TransformationMatrix(ply.Theta);
    sigmaTmp = T*[sigmaX; sigmaY; tauXY];
    sigma1 = sigmaTmp(1);
    sigma2 = sigmaTmp(2);
    tau12 = sigmaTmp(3);
    isFailure = ~(...
        (sigma1 > -ply.Material.sigmaC1ult) && (sigma1 < ply.Material.sigmaT1ult) &&...
        (sigma2 > -ply.Material.sigmaC2ult) && (sigma2 < ply.Material.sigmaT2ult) &&...
        (tau12 > -ply.Material.tau12ult) && (tau12 < ply.Material.tau12ult));
end

function isFailure = maxStrain(sigmaX, sigmaY, tauXY, ply)
    T = TransformationMatrix(ply.Theta);
    
    xiC1ult = ply.Material.sigmaC1ult/ply.Material.E1;
    xiT1ult = ply.Material.sigmaT1ult/ply.Material.E1;
    xiC2ult = ply.Material.sigmaC2ult/ply.Material.E2;
    xiT2ult = ply.Material.sigmaT2ult/ply.Material.E2;
    gamma12ult = ply.Material.tau12ult/ply.Material.G12;
    
    sigmaTmp = T*[sigmaX; sigmaY; tauXY];
    sigma1 = sigmaTmp(1);
    sigma2 = sigmaTmp(2);
    tau12 = sigmaTmp(3);
    
    xi1 = sigma1/ply.Material.E1 - ply.Material.v12*sigma2/ply.Material.E2;
    xi2 = sigma2/ply.Material.E2 - ply.Material.v12*sigma1/ply.Material.E1;
    gamma12 = tau12 / ply.Material.G12;
    
    isFailure = ~(...
        (xi1 > -xiC1ult) && (xi2 < xiT1ult) &&...
        (xi2 > -xiC2ult) && (xi2 < xiT2ult) &&...
        (gamma12 > -gamma12ult) && (gamma12 < gamma12ult));
end

function isFailure = TsaiHill(sigmaX, sigmaY, tauXY, ply)
    T = TransformationMatrix(ply.Theta);
    sigmaTmp = T*[sigmaX; sigmaY; tauXY];
    sigma1 = sigmaTmp(1);
    sigma2 = sigmaTmp(2);
    tau12 = sigmaTmp(3);
    isFailure = ~((...
        (sigma1 / ply.Material.sigmaT1ult) * (sigma1 / ply.Material.sigmaT1ult) -...
        (sigma1 / ply.Material.sigmaT1ult) * (sigma2 / ply.Material.sigmaT1ult) +...
        (sigma2 / ply.Material.sigmaT1ult) * (sigma2 / ply.Material.sigmaT2ult) +...
        (tau12 / ply.Material.tau12ult) * (tau12 / ply.Material.tau12ult)) < 1);
end

function isFailure = TsaiWu(sigmaX, sigmaY, tauXY, ply)
    T = TransformationMatrix(ply.Theta);
    sigmaTmp = T*[sigmaX; sigmaY; tauXY];
    sigma1 = sigmaTmp(1);
    sigma2 = sigmaTmp(2);
    tau12 = sigmaTmp(3);
    
    H1 = 1/ply.Material.sigmaT1ult - 1/ply.Material.sigmaC1ult;
    H2 = 1/ply.Material.sigmaT2ult - 1/ply.Material.sigmaC2ult;
    H6 = 0;
    H11 = 1 / (ply.Material.sigmaT1ult*ply.Material.sigmaC1ult);
    H22 = 1 / (ply.Material.sigmaT2ult*ply.Material.sigmaC2ult);
    H66 = 1 / (ply.Material.tau12ult*ply.Material.tau12ult);
    
    % Mises-Hencky criterion ([1], pg.157)
    H12 = - 0.5*sqrt(1/(ply.Material.sigmaC1ult*ply.Material.sigmaT1ult*ply.Material.sigmaC2ult*ply.Material.sigmaT2ult));
    
    isFailure = ~((...
        H1*sigma1 + H2*sigma2 + H6*tau12 + ...
        H11*sigma1*sigma1 + H22*sigma2*sigma2 + H66*tau12*tau12 +...
        2*H12*sigma1*sigma2) < 1);
end