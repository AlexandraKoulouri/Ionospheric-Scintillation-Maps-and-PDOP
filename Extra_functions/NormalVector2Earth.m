function v = NormalVector2Earth(x0,y0,z0)
%estimate a unit vector normal to Earth at point (x0,y0,z0);

%F(x,y,z) = (x/a)^2 + (y/a)^2 + (z/b)^2 - 1
S = wgs84Ellipsoid('kilometer'); %Earth - Spheroid
a = S.SemimajorAxis;
b = S.SemiminorAxis;

%estimate gradient of F(x,y,z)
Fx = 2*x0/a; Fy=2*y0/a ; Fz = 2*z0/b;

v = [Fx;Fy;Fz;];
v = v./norm(v);
