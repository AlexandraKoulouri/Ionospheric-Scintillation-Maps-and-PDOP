function [P_intersect] = lineIntersect3D(PA,PB)
% Find intersection point of lines in 3D space, in the least squares sense.
% PA :          Nx3-matrix containing starting point of N lines
% PB :          Nx3-matrix containing end point of N lines
% P_Intersect : Best intersection point of the N lines, in least squares sense.
% distances   : Distances from intersection point to the input lines
% Anders Eikenes, 2012
if size(PB,1)==1
    [x,y,z] = Locate_Sat(PB,PA,200000);
    P_intersect = [x y z];
else
    Si = PB - PA; %N lines described as vectors
    ni = Si ./ (sqrt(sum(Si.^2,2))*ones(1,3)); %Normalize vectors
    nx = ni(:,1); ny = ni(:,2); nz = ni(:,3);
    SXX = sum(nx.^2-1);
    SYY = sum(ny.^2-1);
    SZZ = sum(nz.^2-1);
    SXY = sum(nx.*ny);
    SXZ = sum(nx.*nz);
    SYZ = sum(ny.*nz);
    S = [SXX SXY SXZ;SXY SYY SYZ;SXZ SYZ SZZ];
    CX  = sum(PA(:,1).*(nx.^2-1) + PA(:,2).*(nx.*ny)  + PA(:,3).*(nx.*nz));
    CY  = sum(PA(:,1).*(nx.*ny)  + PA(:,2).*(ny.^2-1) + PA(:,3).*(ny.*nz));
    CZ  = sum(PA(:,1).*(nx.*nz)  + PA(:,2).*(ny.*nz)  + PA(:,3).*(nz.^2-1));
    C   = [CX;CY;CZ];
    
    if rcond(S)>eps*10
        P_intersect = (S\C)';
    else
         [x,y,z] = Locate_Sat(PB(1,:),PA(1,:),200000);
    P_intersect = [x y z]; 
        
    end
    
    
end

    

% if nargout>1
%     N = size(PA,1);
%     distances=zeros(N,1);
%     for i=1:N %This is faster:
%         ui=(P_intersect-PA(i,:))*Si(i,:)'/(Si(i,:)*Si(i,:)');
%         distances(i)=norm(P_intersect-PA(i,:)-ui*Si(i,:));
%     end
%     %for i=1:N %http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html:
%     %    distances(i) = norm(cross(P_intersect-PA(i,:),P_intersect-PB(i,:))) / norm(Si(i,:));
%     %end
% end
% end


function [x,y,z,lat,lon] = Locate_Sat(IPP,GroundLoc,h)


%inputs
%IPP - locations of IPP in (x,y,z)
%GroundLon - Location of the received on the ground
%Elevation of the satellite h in km

%Return 
%(x,y,z) coordinates of the satellite
%(lat,lon) of the satellite


% here we have to compute the intersection between the ray and the spheroid
% (Earth)

%this code was created by A. Koulouri, 11.6.2018

%Ray equation
%t --> (x0,y0,z0) + t*(u,v,w)
x0 = GroundLoc(1); y0 = GroundLoc(2); z0 = GroundLoc(3);
u = IPP(1)-GroundLoc(1); v= IPP(2)-GroundLoc(2); w = IPP(3)-GroundLoc(3);


%(x/a)^2 + (y/a)^2 + (z/b)^2 - 1
S = CreateIonosphericSpheroid (h); %satellite trajectory spheroid
a = S.SemimajorAxis^2;
b = S.SemiminorAxis^2;

c1 = b*(x0^2+y0^2)+a*z0^2-a*b;
c2 = 2*b*(u*x0+v*y0)+2*a*w*z0;
c3 = b*(u^2+v^2)+a*w^2;

t1= (-c2+sqrt(c2^2-4*c3*c1))/(2*c3);
t2= (-c2-sqrt(c2^2-4*c3*c1))/(2*c3);
t = max(t1,t2);

% A = 4 *(b^2*(u*x0 + v*y0) + a^2*w*z0)^2 - 4*(b^2*(u^2 + v^2) + a^2*w^2)*(b^2*(-a^2 + x0^2 + y0^2) + a^2*z0^2);
% t = -(1/(b^2*(u^2 + v^2) +  a^2*w^2)) * (b^2*(u*x0+ v*y0) + a^2* w*z0 + 1/2*sqrt(A));

x = x0+t*u;
y = y0+t*v;
z = z0+t*w;


[lat,lon,~] = ecef2geodetic(S,x,y,z);
