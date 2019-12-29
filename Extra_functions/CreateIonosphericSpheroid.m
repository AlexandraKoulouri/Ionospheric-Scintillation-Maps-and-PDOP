function S = CreateIonosphericSpheroid (Altitude)

S =  oblateSpheroid;

EarthS =  wgs84Ellipsoid('kilometer');

S.SemimajorAxis = EarthS.SemimajorAxis+Altitude;
S.SemiminorAxis = EarthS.SemiminorAxis+Altitude;
%S.MeanRadius = (S.SemimajorAxis+S.SemimajorAxis)/2;


%S. Flattening = 1- S.SemiminorAxis/S.SemimajorAxis;
S.Eccentricity = sqrt(2*S.Flattening-S.Flattening^2);
S.InverseFlattening = 1/S.Flattening;
%S.LengthUnit = EarthS.LengthUnit;