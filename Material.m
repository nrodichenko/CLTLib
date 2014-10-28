classdef Material
    properties
        Name
        
        % Source: Table 2.1, p.106, [1]
        % Material properties
Vf          % [unitless] Fiber volume fraction
E1          % [MPa] Longitudinal elastic modulus
E2          % [MPa] Transverse elastic modulus
G12         % [MPa] Shear modulus
v12         % [unitless] Major Poisson ratio

sigmaC1ult  % [MPa] Ultimate longitudinal compressive strength
sigmaT1ult  % [MPa] Ultimate longitudinal tensile strength
sigmaC2ult  % [MPa] Ultimate transverse compressive strength
sigmaT2ult  % [MPa] Ultimate transverse tensile strength
tau12ult    % [MPa] Ultimate in-plane shear strength
        
        % TODO: add the following parameters
        % Longitudinal coefficient of thermal expansion
        % Transverse coefficient of thermal expansion
        % Longitudinal coefficient of moisture expansion
        % Transverse coefficient of moisture expansion
    end
    methods
        % Main constructor
        function m = Material(nName, nVf, nE1, nE2, nG12, nv12, nsigmaC1ult, nsigmaT1ult, nsigmaC2ult, nsigmaT2ult, ntau12ult)
            m.Name = nName;
            m.Vf = nVf;
            m.E1 = nE1;
            m.E2 = nE2;
            m.G12 = nG12;
            m.v12 = nv12;
            m.sigmaC1ult = nsigmaC1ult;
            m.sigmaT1ult = nsigmaT1ult;
            m.sigmaC2ult = nsigmaC2ult;
            m.sigmaT2ult = nsigmaT2ult;
            m.tau12ult = ntau12ult;
        end
        

        
        
    end
    % Static class methods to contain pre-set materials
    methods(Static)
        % Create Glass/Epoxy material
        % Source: Table 2.1, p.106, [1]
        function m = GlassEpoxy()
            m = Material('Glass/Epoxy', 0.45, 38.6e3, 8.27e3, 4.14e3, 0.26, 610, 1062, 118, 31, 72);
        end
        
        % Create Graphite/Epoxy material
        % Source: Table 2.1, p.106, [1]
        function m = GraphiteEpoxy()
            m = Material('Graphite/Epoxy', 0.70, 181e3, 10.30e3, 7.17e3, 0.28, 1500, 1500, 246, 40, 68);
        end
        
        % Create Boron/Epoxy material
        % Source: Table 2.1, p.106, [1]
        function m = BoronEpoxy()
            m = Material('Boron/Epoxy', 0.50, 204e3, 18.50e3, 5.59e3, 0.23, 2500, 1260, 202, 61, 67);
        end
    end
end

