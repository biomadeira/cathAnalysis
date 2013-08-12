from numpy import matrix, sin, cos
Z = lambda theta : matrix([[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0,0,1]])
Y = lambda phi   : matrix([[cos(phi),0, sin(phi)],[0,1,0],[-sin(phi), 0, cos(phi)]])
X = lambda alpha : matrix([[1,0,0],[0, cos(alpha), -sin(alpha)],[0, sin(alpha), cos(alpha)]])

#euler rotationMatrices Convention
XYZ = lambda alpha, phi, theta : X(alpha)*Y(phi)*Z(theta)
ZXY = lambda alpha, phi, theta : Z(theta)*X(alpha)*Y(phi)
YXZ = lambda alpha, phi, theta : Y(theta)*X(alpha)*Z(phi)

ZXZ = lambda theta, alpha, theta2 : Z(theta)*X(alpha)*Z(theta2)
ZYZ = lambda theta, phi, theta2   : Z(theta)*Y(phi)*Z(theta2)
ZYX = lambda theta, phi, alpha    : Z(theta)*Y(phi)*X(alpha)

transform = lambda atom_coords, center, shift, rotationMatrix: rotationMatrix*(atom_coords - center)+shift