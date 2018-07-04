function [ Bi,Gi,A_tetai,B_bari ] = Matrices( DNX,DNY,i,U )
    Gi=zeros(4,2);
    Bi=zeros(3,2);
    A_tetai=zeros(3,4);
    B_bari=zeros(3,2);
    Gi=[DNX(i,1),0;0,DNX(i,1);DNY(i,1),0;0,DNY(i,1)];
    Bi=[DNX(i,1),0;0,DNY(i,1);DNY(i,1),DNX(i,1)];
    for j=1:4
        A_tetai=[DNX(j,1)*U(2*j-1,1),DNX(j,1)*U(2*j,1),0,0;0,0,DNY(j,1)*U(2*j-1,1),DNY(j,1)*U(2*j,1);DNY(j,1)*U(2*j-1,1),DNY(j,1)*U(2*j,1),DNX(j,1)*U(2*j-1,1),DNX(j,1)*U(2*j,1)]+A_tetai;
    end
    B_bari=Bi+A_tetai*Gi;

end

