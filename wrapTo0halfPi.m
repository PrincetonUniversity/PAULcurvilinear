function output=wrapTo0halfPi(angleInRad)
% the input needs to be an angle in [0, pi]
angleInRad=wrapTo0Pi(angleInRad);
coeff=((0.5.*pi-angleInRad)>=0);
output=(1-coeff).*pi+(coeff.*2-1).*angleInRad;

function output=wrapTo0Pi(angleInRad)

temp = wrapToPi(angleInRad);
output = temp+pi*(temp<0);
