function output=wrapTo0halfPi(angleInRad)
% the input needs to be an angle in [0, pi]
angleInRad=wrapTo0Pi(angleInRad);
coeff=sign(0.5.*pi-angleInRad);
output=((1-coeff)/2).*pi+coeff.*angleInRad;

function output=wrapTo0Pi(angleInRad)

temp = wrapToPi(angleInRad);
output = temp+pi*((1-sign(temp))/2);