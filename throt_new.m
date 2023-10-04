function [M] = throt_new(alpha,phi)

ca = cos(pi/180*alpha);	% cosine of tip alpha
sa = sin(pi/180*alpha);	% sine of tip
cp = cos(pi/180*phi  ); % cosine of phi
sp = sin(pi/180*phi  ); % sine of phi


M = [cp*cp+sp*sp*ca cp*sp*(1-ca) -sp*sa;
     cp*sp-sp*cp*ca sp*sp+cp*cp*ca cp*sa;
     sa*sp -sa*cp ca];
end