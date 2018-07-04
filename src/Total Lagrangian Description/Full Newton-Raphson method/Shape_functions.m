function [ N,DNX,DNY,jacobian] = Shape_functions( xi,yi,Element_coordinate )

N(1,1)=0.25*(1-xi)*(1+yi);
N(2,1)=0.25*(1-xi)*(1-yi);
N(3,1)=0.25*(1+xi)*(1-yi);
N(4,1)=0.25*(1+xi)*(1+yi);

DDNX(1,1)=-0.25*(1+yi);
DDNX(2,1)=-0.25*(1-yi);
DDNX(3,1)=0.25*(1-yi);
DDNX(4,1)=0.25*(1+yi);

DDNY(1,1)=0.25*(1-xi);
DDNY(2,1)=-0.25*(1-xi);
DDNY(3,1)=-0.25*(1+xi);
DDNY(4,1)=0.25*(1+xi);

jacobian=zeros(2,2);
for i=1:4;
    x=Element_coordinate(i,2);
    y=Element_coordinate(i,3);
    jacobian(1,1)=jacobian(1,1)+x*DDNX(i,1);
    jacobian(1,2)=jacobian(1,2)+y*DDNX(i,1);
    jacobian(2,1)=jacobian(2,1)+x*DDNY(i,1);
    jacobian(2,2)=jacobian(2,2)+y*DDNY(i,1);
end

for i=1:4;
    DNX(i,1)=(jacobian(2,2)*DDNX(i,1)-jacobian(1,2)*DDNY(i,1))/det(jacobian);
    DNY(i,1)=(jacobian(1,1)*DDNY(i,1)-jacobian(2,1)*DDNX(i,1))/det(jacobian);
end

end