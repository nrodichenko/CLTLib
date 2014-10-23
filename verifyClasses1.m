function out = verifyClasses1()

stresses = zeros(4, 3);
stresses(1, :) = [100.0,	0.0,	0.0];
stresses(2, :) = [0.0,	50.0,	0.0];		
stresses(3, :) = [0.0,	0.0,	100.0];
stresses(4, :) = [100.0,	0.0,	-100.0];
stresses = stresses * 1e6;
NS = size(stresses, 1);

materialSystems = zeros(1, 4);
materialSystems(1, :) = [181.0,	10.3,	7.17,	0.28];

thetas = [0, pi/8, pi/4, pi/3];
NT = size(thetas, 2);

out = cell(2, 1);

for k=1:1
	t = cell(NS, NT*3);
	tt = cell(NS, NT*3);
	for i=1:NS
		for j=1:NT
			t1 = getStrainsFromStress(stresses(i, 1), stresses(i, 2), stresses(i, 3), materialSystems(k, 1), materialSystems(k, 2), materialSystems(k, 3), materialSystems(k, 4), thetas(j));
			t{i, (j-1)*3 + 1} = t1(1);
			t{i, (j-1)*3 + 2} = t1(2);
			t{i, (j-1)*3 + 3} = t1(3);
            
            ply = Ply(Material.GraphiteEpoxy(), thetas(j));
			t2 = ply.GetStrainsFromStress(stresses(i, :));
			tt{i, (j-1)*3 + 1} = t2(1);
			tt{i, (j-1)*3 + 2} = t2(2);
			tt{i, (j-1)*3 + 3} = t2(3);
		end
	end
	out{k} = t;
	out{k+1} = tt;
end



end


function out = getStrainsFromStress(sigmaX, sigmaY, tauXY, E1, E2, G12, v12, theta)

Qd = Qdash(E1, E2, G12, v12, theta);

out = inv(Qd) * [sigmaX; sigmaY; tauXY];

eps1 = out(1);
eps2 = out(2);
gamma12 = out(3);

end

function out = Qdash(E1, E2, G12, v12, theta)

v21 = E2 * v12 / E1;
Q11 = E1 / (1 - v12 * v21);
Q22 = E2 / (1 - v12 * v21);
Q12 = v12 * E2 / (1 - v12 * v21);
Q66 = G12;

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