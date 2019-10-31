clear
clc
%INTERVAL Ci AND Cij:
Liminf=0;
LimSup=2;
d=0.05;
e=1./d;
k=(LimSup./d).^2;  %Number of rows%
%MATRIX Cij%
aa1=3;
aa2=4;
s=linspace(Liminf,LimSup,LimSup./d);
H=s;
for i=1:(LimSup/d)-1
    C12=[s,H];
    H=C12;
end
count16=1;
for i=1:k
       if H(1,i)<LimSup
       C21(1,i)=H(1,count16);
       else
           C21(1,i)=H(1,count16);
           count16=count16+1;
          end
end
for i=1:k
    C1(i,1)=aa1;
end
for i=1:k
    C2(i,1)=aa2;
end
Coefacpp=[C1,C2,H',C21'];
[ee rr]=size(Coefacpp);

for i=1:ee
    Coefacp(i,1)=prod(Coefacpp(i,[1 3 4]));%(only C2>=0)% 
end
Coefacz=[Coefacpp,Coefacp];
count70=1;
for i=1:ee
    if Coefacz(i,5)>0
        Coefac(count70,:)=Coefacz(i,:);
        count70=count70+1;
    end
end

%_______________________________________________________________________________________________________________%
%COEFFICIENT GAMMA OF GENERAL POLYNOMIAL G(Ci,Cij):
Ga9=((Coefac(:,1).^3).*Coefac(:,2))./(Coefac(:,3).^3);
Ga8=((-3).*(Coefac(:,1).^3).*Coefac(:,2))./(Coefac(:,3).^3);
Ga7=(3.*(Coefac(:,1).^2).*Coefac(:,2).*(1+Coefac(:,1)))./(Coefac(:,3).^3);
Ga6=((-1).*(Coefac(:,1).^2).*Coefac(:,2).*(-Coefac(:,3)+Coefac(:,1)+9))./(Coefac(:,3).^3);
Ga5=((-1).*Coefac(:,1).*Coefac(:,2).*(Coefac(:,1).*(-9+2.*Coefac(:,3))-3))./(Coefac(:,3).^3);
Ga4=(Coefac(:,1).*Coefac(:,2).*(Coefac(:,1).*Coefac(:,3)-3.*Coefac(:,1)+2.*Coefac(:,3)-9))./(Coefac(:,3).^3);
Ga3=(-4.*Coefac(:,1).*Coefac(:,2).*Coefac(:,3)+Coefac(:,1).*Coefac(:,3).^2+9.*Coefac(:,1).*Coefac(:,2)+Coefac(:,2))./(Coefac(:,3).^3);
Ga2=(2.*Coefac(:,1).*Coefac(:,2).*Coefac(:,3)-Coefac(:,1).*Coefac(:,3).^2-3.*Coefac(:,1).*Coefac(:,2)+Coefac(:,2).*Coefac(:,3)-3.*Coefac(:,2))./(Coefac(:,3).^3);
Ga1=(-Coefac(:,4).*Coefac(:,3).^3-2.*Coefac(:,3).*Coefac(:,2)+Coefac(:,3).^2+3.*Coefac(:,2))./(Coefac(:,3).^3);
Ga0=(Coefac(:,3)-1).*(Coefac(:,3).^2+Coefac(:,2))./(Coefac(:,3).^3);
format long
Mgamma=[Ga9,Ga8,Ga7,Ga6,Ga5,Ga4,Ga3,Ga2,Ga1,Ga0];
%___________________________________________________________________________________________________________________%
%MATRIZ DE CEROS DEL POLONOMIO GENERAL EN COMPETICION POR INT:
[cc vv]=size(Mgamma);
for j= 1:cc
Pn=Mgamma(j,:);
RootPn=roots(Pn);
S(j,:)=RootPn;
end
format long
S;

% % % %________________________________________________________________________________________________________________%
% %POSITIVE ROOTS OF POLINOMIAL
[z,o]=size(S);
count=0;
for i=1:z
    count25=0;
        for n=1:o
           if imag(S(i,n))==0
             if real(S(i,n))>0
              %if real(S(i,n))<1   
               count25=count25+1;
               count=count+1;
               CS(i,1)=count25; 
             LPRTri(i,n)=S(i,n);          
             end
             end
           end
        end
          % end
format short
% % % %**********************************************Stability*****************************************************%
     LPRTricount=[LPRTri,Coefac,CS];
     [wq,en]=size(LPRTricount);
     if aa1.*aa2==0
   LTRCcol=[LPRTricount(:,1),LPRTricount(:,[4:7 9]);LPRTricount(:,2),LPRTricount(:,[4:7 9]);LPRTricount(:,3),LPRTricount(:,[4:7 9])]; 
         else
LTRCcol=[LPRTricount(:,1),LPRTricount(:,[10:13 15]);LPRTricount(:,2),LPRTricount(:,[10:13 15]);LPRTricount(:,3),LPRTricount(:,[10:13 15]);LPRTricount(:,4),LPRTricount(:,[10:13 15]);LPRTricount(:,5),LPRTricount(:,[10:13 15]);LPRTricount(:,6),LPRTricount(:,[10:13 15]);LPRTricount(:,7),LPRTricount(:,[10:13 15]);LPRTricount(:,8),LPRTricount(:,[10:13 15]);LPRTricount(:,9),LPRTricount(:,[10:13 15])];     
     end
[cw,mp]=size(LTRCcol); 
 count28=0;       
for i=1:cw
    if LTRCcol(i,1)>0
        count28=count28+1;
        LTRCC(count28,:)=LTRCcol(i,:);
    end
end
% %Function Nullclines (U1,U2)%
% %________________________________________________________________________________________________________________%
[ap,pe]=size(LTRCC);
for i=1:ap
    POR1(i,1)=LTRCC(i,1);
    POR2(i,1)=(1./LTRCC(i,4)).*(-LTRCC(i,2).*POR1(i,1).^3+LTRCC(i,2).*POR1(i,1).^2-POR1(i,1)+1); 
end
LTRCU=[POR2,LTRCC];
count18=0;
for i=1:ap
    if LTRCU(i,1)>0
            count18=count18+1;
        PORRC(count18,:)=LTRCU(i,:);
        end
end
Uscoef=PORRC;
[er,ty]=size(Uscoef);
% %-------Jacobian Matrix-----------%%%%%%
syms u1 u2
for i=1:er  
cc1=Uscoef(i,3);
cc2=Uscoef(i,4);
cc12=Uscoef(i,5);
cc21=Uscoef(i,6);
F=[u1-u1.^2-((cc12.*u1.*u2)./(1+cc1.*u1.^2)),u2-u2.^2-((cc21.*u1.*u2)./(1+cc2.*u2.^2))];
v=[u1 u2];
J=jacobian(F,v);
PP=[Uscoef(i,2) Uscoef(i,1)];
JN=double(subs(J,v,PP));
Lanc(i,:)=eig(JN)';
end
Lan=[Lanc,Uscoef];
% % %FORCE AND ABILITY TO LIMIT SPECIES TO STABILITY%
% % %________________________________________________________________________________________________________________%
[qw,as]=size(Lan);
count3=0;
for i=1:qw
        if Lan(i,1)<0
            if Lan(i,2)<0
                count3=count3+1;
                Lanneg(count3,:)=Lan(i,:);
            end
        end     
end
[m,n]=size(Lanneg);
% %graphics Cij%
% %________________________________________________________________________________________________________________%
%%Data matrix of coexistence 0<Cij<1%%
count30=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)<1
        if Lanneg(i,9)==1
        count30=count30+1;
        coe1(count30,1)=Lanneg(i,7);
        coe2(count30,:)=Lanneg(i,[8 9]);
        coe011=[coe1,coe2];
        end
    end
end
end

count31=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)<1
        if Lanneg(i,9)==3
        count31=count31+1;
        coe13(count31,1)=Lanneg(i,7);
        coe23(count31,:)=Lanneg(i,[8 9]);
        coe013=[coe13,coe23];
        end
    end
end
end
count32=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)<1
        if Lanneg(i,9)==5
        count32=count32+1;
        coe15(count32,1)=Lanneg(i,7);
        coe25(count32,:)=Lanneg(i,[8 9]);
        coe015=[coe15,coe25];
        end
    end
end
end
%%________________________________________________________________________________________%
%%Data matrix of Exclusion Cij>1%%
count33=0;
for i=1:m
if Lanneg(i,7)>1
    if Lanneg(i,8)>1
        if Lanneg(i,9)==3
        count33=count33+1;
        exclu1(count33,1)=Lanneg(i,7);
        exclu2(count33,:)=Lanneg(i,[8 9]);
        exclu=[exclu1,exclu2];
        end
    end
end
end
%%________________________________________________________________________________________%%
%%Data matrix of SP1 C12<1 and C21>1%%
count34=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)>1
        if Lanneg(i,9)==2
        count34=count34+1;
        sp112(count34,1)=Lanneg(i,7);
        sp122(count34,:)=Lanneg(i,[8 9]);
        sp12=[sp112,sp122];
        end
    end
end
end
count35=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)>1
        if Lanneg(i,9)==4
        count35=count35+1;
        sp114(count35,1)=Lanneg(i,7);
        sp124(count35,:)=Lanneg(i,[8 9]);
        sp14=[sp114,sp124];
        end
    end
end
end
%%_____________________________________________________________________________________%%
%%Data matrix of SP1 C12>1 and C21<1%%
count36=0;
for i=1:m
if Lanneg(i,7)>1
    if Lanneg(i,8)<1
        if Lanneg(i,9)==2
        count36=count36+1;
        sp212(count36,1)=Lanneg(i,7);
        sp222(count36,:)=Lanneg(i,[8 9]);
        sp22=[sp212,sp222];
        end
    end
end
end
count37=0;
for i=1:m
if Lanneg(i,7)>1
    if Lanneg(i,8)<1
        if Lanneg(i,9)==4
        count37=count37+1;
        sp214(count37,1)=Lanneg(i,7);
        sp224(count37,:)=Lanneg(i,[8 9]);
        sp24=[sp214,sp224];
        end
    end
end
end
%%__________________________________________________________________________________________________%%
%%FIGURE%%
figure
if count30>0
scatter(coe011(:,1),coe011(:,2),'.','g')
end
hold on
if count31>0
scatter(coe013(:,1),coe013(:,2),'.','y')
end
if count32>0
scatter(coe015(:,1),coe015(:,2),'.','k')
end
if count34>0
scatter(sp12(:,1),sp12(:,2),'.','r')
end
if count35>0
scatter(sp14(:,1),sp14(:,2),'.','m')
end

if count33>0
scatter(exclu(:,1),exclu(:,2),'+','k')
end
if count36>0
scatter(sp22(:,1),sp22(:,2),'x','b')
end
if count37>0
scatter(sp24(:,1),sp24(:,2),'x','c')
end

ylim([0 8])
xlim([0 8])
xlabel('c_{12}'),ylabel('c_{21}')
hold off
grid on
grid minor
