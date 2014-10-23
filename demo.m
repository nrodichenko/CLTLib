% Define a ply with pre-set material
ply = Ply(Material.BoronEpoxy(), pi/3);

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
ply.DrawFailureEnvelope(0, 'Tsai-Wu'); 



% Draw failure envelope for tauXY=0 with MaximumStrength criterion
ply.DrawFailureEnvelope(0, 'MaximumStrength'); 