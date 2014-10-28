classdef Ply
    properties
        Material    % Ply material. Instance of Material class.
        Theta       % [radians]Ply angle
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
        % Possible criterions: MaximumStress, MaximumStrain, Tsai-Hill, Tsai-Wu
        function out = CheckFailure(obj, StressVector, FailureCriterion)
            sigmaX = StressVector(1);
            sigmaY = StressVector(2);
            tau = StressVector(3);
            
            switch FailureCriterion
                case 'MaximumStress'
                    fcn = MaxStressFcn(sigmaX, sigmaY, tau, obj);
                case 'MaximumStrain'
                    fcn = MaxStrainFcn(sigmaX, sigmaY, tau, obj);
                case 'Tsai-Hill'
                    fcn = TsaiHillFcn(sigmaX, sigmaY, tau, obj);
                case 'Tsai-Wu'
                    fcn = TsaiWuFcn(sigmaX, sigmaY, tau, obj);
                otherwise
                    error([FailureCriterion ' is not a valid failure criterion. Possible criterions: MaximumStress, MaximumStrain, Tsai-Hill, Tsai-Wu'])
            end 
            
            out = ~(fcn < 0);
        end        
        
        
        % Draw failure envelope for given tauXY
        % and failure criterion
        % Possible criterions: MaximumStress, MaximumStrain, Tsai-Hill, Tsai-Wu
        function DrawFailureEnvelope(obj, tau, FailureCriterion)
            figure;
            if ~iscell(FailureCriterion)
               tmp = cell(1,1);
               tmp{1} = FailureCriterion;
               FailureCriterion = tmp;
            end
            
            hold all;
            colormap = lines(length(FailureCriterion));
            for i = 1:length(FailureCriterion)
               fc = FailureCriterion{i};
               
                switch fc
                    case 'MaximumStress'
                        fcn = @(sigmaX, sigmaY)Ply.MaxStressFcn(sigmaX, sigmaY, tau, obj);
                    case 'MaximumStrain'
                        fcn = @(sigmaX, sigmaY)Ply.MaxStrainFcn(sigmaX, sigmaY, tau, obj);
                    case 'Tsai-Hill'
                        fcn = @(sigmaX, sigmaY)Ply.TsaiHillFcn(sigmaX, sigmaY, tau, obj);
                    case 'Tsai-Wu'
                        fcn = @(sigmaX, sigmaY)Ply.TsaiWuFcn(sigmaX, sigmaY, tau, obj);

                    otherwise
                        error([FailureCriterion ' is not a valid failure criterion. Possible criterions: MaximumStress, MaximumStrain, Tsai-Hill, Tsai-Wu'])
                end

                h(i) = ezplot(fcn, [-1e3,1e3,-1e3,1e3]);
                set(h(i), 'Color', colormap(i, :));
            end
            
            set(gca, 'XMinorTick','on','FontWeight','bold',...
                'FontSize',14,'YGrid','on','XGrid','on');
            ylabel('\sigma_y [MPa]');
            xlabel('\sigma_x [MPa]');
            title({...
                    ['Failure envelopes for single-ply ' obj.Material.Name ', \theta = ' num2str(obj.Theta*180.0/pi, '%.1f') char(176)]});
                
            legend(FailureCriterion);
            hold off;            
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
        
        function out = Q(obj)
            v21 = obj.Material.E2 * obj.Material.v12 / obj.Material.E1;
            Q11 = obj.Material.E1 / (1 - obj.Material.v12 * v21);
            Q22 = obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q12 = obj.Material.v12 * obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q66 = obj.Material.G12;
            out = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
        end
        
        % Seems to be incorrect. But why?
        function out = QBarAlt(obj)
            v21 = obj.Material.E2 * obj.Material.v12 / obj.Material.E1;
            Q11 = obj.Material.E1 / (1 - obj.Material.v12 * v21);
            Q22 = obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q12 = obj.Material.v12 * obj.Material.E2 / (1 - obj.Material.v12 * v21);
            Q66 = obj.Material.G12;
            Qmod = [Q11 Q12 0; Q12 Q22 0; 0 0 2*Q66];
            
            T = obj.TransformationMatrix();
            out = inv(T)*Qmod*T;
        end
        
        
        % Calculates strains from provided stresses and ply properties
        function strains = GetStrainsFromStress(obj, StressVector)
            qb = obj.QBar();
            strains = inv(qb)*reshape(StressVector, [3,1]); % The same as inv(qb) * StressVector
        end
    end
    
    methods(Static)

        % Helper methods

        function out = MaxStressFcn(sigmaXvec, sigmaYvec, tauXY, ply)
            out = zeros(length(sigmaXvec), 1);

            T = ply.TransformationMatrix();

            for i=1:length(sigmaXvec)
                sigmaX = sigmaXvec(i);
                sigmaY = sigmaYvec(i);
                sigmaTmp = T*[sigmaX; sigmaY; tauXY];
                sigma1 = sigmaTmp(1);
                sigma2 = sigmaTmp(2);
                tau12 = sigmaTmp(3);

                % below we construct a continuos function that has zero value on
                % the failure envelope boundary, is less then zero outside failure
                % envelope, and is greater then zero inside the failure envelope

                % we need such a function to be able to feed the constraints into
                % the Matlab's ezplot routine

                % first, to construct a continuous function we take the minimum
                % of squared distances to each boundary
                out(i) = ...
                    min((sigma1 + ply.Material.sigmaC1ult)^2, ...
                    min((sigma1 - ply.Material.sigmaT1ult)^2, ...
                    min((sigma2 + ply.Material.sigmaC2ult)^2, ... 
                    min((sigma2 - ply.Material.sigmaT2ult)^2, ...
                    min((tau12 + ply.Material.tau12ult)^2, ...
                        (tau12 - ply.Material.tau12ult)^2 ...
                    )))));

                % isFailure == true outside failure envelope
                isFailure = ~(...
                    (sigma1 > -ply.Material.sigmaC1ult) && (sigma1 < ply.Material.sigmaT1ult) &&...
                    (sigma2 > -ply.Material.sigmaC2ult) && (sigma2 < ply.Material.sigmaT2ult) &&...
                    (tau12 > -ply.Material.tau12ult) && (tau12 < ply.Material.tau12ult));

                % constructed function is positive outside failure envelope, and
                % negative inside the failure envelope
                if isFailure
                    sign = +1;
                else
                    sign = -1;
                end

                out(i) = out(i) * sign;
            end
        end

        function out = MaxStrainFcn(sigmaXvec, sigmaYvec, tauXY, ply)
            out = zeros(length(sigmaXvec), 1);

            T = ply.TransformationMatrix();

            xiC1ult = ply.Material.sigmaC1ult/ply.Material.E1;
            xiT1ult = ply.Material.sigmaT1ult/ply.Material.E1;
            xiC2ult = ply.Material.sigmaC2ult/ply.Material.E2;
            xiT2ult = ply.Material.sigmaT2ult/ply.Material.E2;
            gamma12ult = ply.Material.tau12ult/ply.Material.G12;


            for i=1:length(sigmaXvec)
                sigmaX = sigmaXvec(i);
                sigmaY = sigmaYvec(i);
                sigmaTmp = T*[sigmaX; sigmaY; tauXY];
                sigma1 = sigmaTmp(1);
                sigma2 = sigmaTmp(2);
                tau12 = sigmaTmp(3);

                xi1 = sigma1/ply.Material.E1 - ply.Material.v12*sigma2/ply.Material.E2;
                xi2 = sigma2/ply.Material.E2 - ply.Material.v12*sigma1/ply.Material.E1;
                gamma12 = tau12 / ply.Material.G12;

                % below we construct a continuos function that has zero value on
                % the failure envelope boundary, is less then zero outside failure
                % envelope, and is greater then zero inside the failure envelope

                % we need such a function to be able to feed the constraints into
                % the Matlab's ezplot routine

                % first, to construct a continuous function we take the minimum
                % of squared distances to each boundary
                out(i) = ...
                    min((xi1 + xiC1ult)^2, ...
                    min((xi1 - xiT1ult)^2, ...
                    min((xi2 + xiC2ult)^2, ...
                    min((xi2 - xiT2ult)^2, ...
                    min((gamma12 - gamma12ult)^2, ...
                        (gamma12 + gamma12ult)^2 ...
                    )))));

                % isFailure == true outside failure envelope
                isFailure = ~(...
                    (xi1 > -xiC1ult) && (xi2 < xiT1ult) &&...
                    (xi2 > -xiC2ult) && (xi2 < xiT2ult) &&...
                    (gamma12 > -gamma12ult) && (gamma12 < gamma12ult));

                % constructed function is positive outside failure envelope, and
                % negative inside the failure envelope
                if isFailure
                    sign = +1;
                else
                    sign = -1;
                end

                out(i) = out(i) * sign;
            end
        end

        function out = TsaiHillFcn(sigmaXvec, sigmaYvec, tauXY, ply)    
            out = zeros(length(sigmaXvec), 1);

            T = ply.TransformationMatrix();
            for i=1:length(sigmaXvec)
                sigmaX = sigmaXvec(i);
                sigmaY = sigmaYvec(i);

                sigmaTmp = T*[sigmaX; sigmaY; tauXY];
                sigma1 = sigmaTmp(1);
                sigma2 = sigmaTmp(2);
                tau12 = sigmaTmp(3);

                out(i) = (...
                (sigma1 / ply.Material.sigmaT1ult) * (sigma1 / ply.Material.sigmaT1ult) -...
                (sigma1 / ply.Material.sigmaT1ult) * (sigma2 / ply.Material.sigmaT1ult) +...
                (sigma2 / ply.Material.sigmaT1ult) * (sigma2 / ply.Material.sigmaT2ult) +...
                (tau12 / ply.Material.tau12ult) * (tau12 / ply.Material.tau12ult)) - 1;
            end
        end

        function out = TsaiWuFcn(sigmaXvec, sigmaYvec, tauXY, ply)
            out = zeros(length(sigmaXvec), 1);

            T = ply.TransformationMatrix();

            H1 = 1/ply.Material.sigmaT1ult - 1/ply.Material.sigmaC1ult;
            H2 = 1/ply.Material.sigmaT2ult - 1/ply.Material.sigmaC2ult;
            H6 = 0;
            H11 = 1 / (ply.Material.sigmaT1ult*ply.Material.sigmaC1ult);
            H22 = 1 / (ply.Material.sigmaT2ult*ply.Material.sigmaC2ult);
            H66 = 1 / (ply.Material.tau12ult*ply.Material.tau12ult);
            % Mises-Hencky criterion ([1], pg.157)
            H12 = - 0.5*sqrt(1/(ply.Material.sigmaC1ult*ply.Material.sigmaT1ult*ply.Material.sigmaC2ult*ply.Material.sigmaT2ult));

            for i=1:length(sigmaXvec)
                sigmaX = sigmaXvec(i);
                sigmaY = sigmaYvec(i);

                sigmaTmp = T*[sigmaX; sigmaY; tauXY];
                sigma1 = sigmaTmp(1);
                sigma2 = sigmaTmp(2);
                tau12 = sigmaTmp(3);

                out(i) = (...
                    H1*sigma1 + H2*sigma2 + H6*tau12 + ...
                    H11*sigma1*sigma1 + H22*sigma2*sigma2 + H66*tau12*tau12 +...
                    2*H12*sigma1*sigma2) - 1;
            end
        end
        
    end
end

