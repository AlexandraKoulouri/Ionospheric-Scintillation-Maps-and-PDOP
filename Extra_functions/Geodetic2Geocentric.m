function  [x,y,z] = Geodetic2Geocentric(lon,lat,Alt,Spher)
%this code converts the geodetic coordinates to geocentric (Cartesian)

%this code was created by A. Koulouri, 26.5.2018
if nargin == 3
 Spher = wgs84Ellipsoid('kilometer'); %Earth spheroid
end
%Convert geodetic to geocentric (ECEF) coordinates
[x,y,z] = geodetic2ecef(Spher,lat,lon,Alt);

