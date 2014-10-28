% Define a ply with pre-set material
mat = Material.GraphiteEpoxy();
ply = Ply(mat, pi/3);

% Define a ply with custom material
customMaterial2 = Material('Graphite/Epoxy-Strong', 0.70, 181, 10.30, 7.17, 0.28, 3000, 3000, 500, 80, 130);
customPly2 = Ply(customMaterial2, pi/3);


% Calculate strains from stresses
stresses = zeros(4, 3);
stresses(1, :) = [100.0,	0.0,	0.0];
stresses(2, :) = [0.0,	50.0,	0.0];		
stresses(3, :) = [0.0,	0.0,	100.0];
stresses(4, :) = [100.0,	0.0,	-100.0];

for i = 1:size(stresses, 1)
    ply.GetStrainsFromStress(stresses(i, :))
end


% Draw failure envelope for tauXY=0 with Tsai-Wu criterion


% Draw failure envelope for tauXY=0 with MaximumStrain criterion
customPly2.DrawFailureEnvelope(0, 'MaximumStrain'); 
ply.DrawFailureEnvelope(0, 'MaximumStrain');  