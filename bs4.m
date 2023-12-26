clear all
clc 
x1=1.5,x2=0.6,x3=-1.5;
jx1=0,jx2=0,jx3=0;
jw01=0,jw02=0,jw03=0; 
jy1=0,jy2=0,jy3=0;  
t=0;
h=0.0001;
alpha=0.8;
A1=[-2.14,-1.14,0;1.114,-2.016,0;1.00,0.15,-2];
A2=[-2.01,1.04,0;0,-0.4,-1.114;1.03,0,-1.01];
B1=[1.25;0;0];
B2=[0;1.35;0];
C1=[1.01,0,0;-1.01,1.01,0;0,0,1.01];%C2=C1;
C2=[1.01,0,0;1.01,0,1.01;0,1.01,0];
C=C1;
lmd1=[4,0,0;0,3,0;0,0,2];
lmd2=[4,0,0;0,3,0;0,0,2];
dlA1=[0.060,0,0;0,0.05,0;0,0,0.02];
dlA2=[0.2215,0,0;0,0.15,0;0,0,0.02];
epsl1=norm(dlA1);
epsl2=norm(dlA2);
gam1=0.05;
gam2=0.12;
hE=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0];
hA1=[-2.14,-1.14,0,0,0,0;1.114,-2.016,0,0,0,0;1,0.15,-2,0,0,0;0,-0,0,-1,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
hA2=[-2.01,1.04,0,0,0,0;0,-0.4,-1.114,0,0,0;1.03,0,-1.01,0,0,0;0,0,0,-1,0,0;0,0,0,0,-1,0;0,0,0,0,0,-1];
hB1=[1.25;0;0;0;0;0];hB2=[0;1.35;0;0;0;0];
hN=[0,0,0;0,0,0;0,0,0;1,0,0;0,1,0;0,0,1];
%hC=[1,0,0,1,0,0;-1,1,0,0,1,0;0,0,1,0,0,1];
%C0=[1,0,0,0,0,0;-1,1,0,0,0,0;0,0,1,0,0,0];
hC=[1.01,0,0,1,0,0;-1.01,1.01,0,0,1,0;0,0,1.01,0,0,1];
C0=[1.01,0,0,0,0,0;-1.01,1.01,0,0,0,0;0,0,1.01,0,0,0];
hL=[0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1]';
hh=inv(hE+hL*hC);
epsl11=norm(hh)*epsl1;
epsl12=norm(hh)*epsl2;
p1=epsl11/(norm(hC)*gam1);
p2=epsl12/(norm(hC)*gam2);
a11=0,a12=0,a13=0,a14=0,a15=0,a16=0,a21=0,a22=0,a23=0,a24=0,a25=0,a26=0,a31=0,a32=0,a33=0,a34=0,a35=0,a36=0;
b11=0,b12=0,b13=0,b14=0,b15=0,b16=0,b21=0,b22=0,b23=0,b24=0,b25=0,b26=0,b31=0,b32=0,b33=0,b34=0,b35=0,b36=0;
jhth1=[a11,a12,a13,a14,a15,a16;a21,a22,a23,a24,a25,a26;a31,a32,a33,a34,a35,a36];
jhth2=[b11,b12,b13,b14,b15,b16;b21,b22,b23,b24,b25,b26;b31,b32,b33,b34,b35,b36];


for n=1:149990
    t(n+1)=t(n)+h;
    u(n)=sin(t(n));
    h1(n)=1-(sin(x1(n)))^2;
    h2(n)=(sin(x1(n)))^2;
    dA1=dlA1*cos(t(n));
    dA2=dlA2*cos(t(n));
    hdA1=[0.06,0,0,0,0,0;0,0.05,0,0,0,0;0,0,0.02,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0]*cos(t(n));
    hdA2=[0.225,0,0,0,0,0;0,0.15,0,0,0,0;0,0,0.02,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0]*cos(t(n));
    jz=[x1(n);x2(n);x3(n)];
    AAA=h1(n)*((A1+dA1)*jz+B1*u(n))+h2(n)*((A2+dA2)*jz+B2*u(n));
    f1(n)=AAA(1);
    f2(n)=AAA(2);
    f3(n)=AAA(3);
    G=(t(n+1)-t(1:n)).^(alpha-1);
    ff1=G.*f1/gamma(alpha);
    ff2=G.*f2/gamma(alpha);
    ff3=G.*f3/gamma(alpha);
    x1(n+1)=x1(1)+h*sum(ff1);
    x2(n+1)=x2(1)+h*sum(ff2);
    x3(n+1)=x3(1)+h*sum(ff3);
    if t(n)>=2
        w1(n)=0.5*(sin(5*(t(n)-5)*pi))^2;
    else
        w1(n)=0;
    end
    if t(n)>=4
        w2(n)=0.2*(cos(2*(t(n)-4)*pi))^2;
    else
        w2(n)=0;
    end
        w3(n)=0.3*((sin(3*pi*t(n)))^2+cos(2*pi*t(n)));
    jzw=[w1(n);w2(n);w3(n)];
    BBB=h1(n)*C1*jz+h2(n)*C2*jz+jzw;
    EEE=h1(n)*(C1-C)*jz+h2(n)*(C2-C)*jz+jzw;
    w01(n)=EEE(1);
    w02(n)=EEE(2);
    w03(n)=EEE(3);
    y1(n)=BBB(1);
    y2(n)=BBB(2);
    y3(n)=BBB(3);
    yy=[y1(n);y2(n);y3(n)];
    jyy=hC*[jx1(n);jx2(n);jx3(n);jw01(n);jw02(n);jw03(n)];
    hA01=hA1-hN*C0;
    hA02=hA2-hN*C0;
    jzjhx0=[jx1(n);jx2(n);jx3(n);jw01(n);jw02(n);jw03(n)];
    yjy=yy-jyy;
    hA11=hh*hA01;
    hA12=hh*hA02;

    setlmis([]);
    P=lmivar(1,[6,1]);
    Q1=lmivar(2,[6,3]);
    Q2=lmivar(2,[6,3]);
    %1st LMI
    lmiterm([1 1 1 P],1,hA11,'s');
    lmiterm([1 1 1 Q1],-1,hC,'s');
    lmiterm([1 1 1 P],2*epsl11,1);
    %2st LMI
    lmiterm([1 1 1 P],1,hA12,'s');
    lmiterm([1 1 1 Q2],-1,hC,'s');
    lmiterm([1 1 1 P],2*epsl12,1);
    lmis=getlmis;
    [tmin,xfeas]=feasp(lmis,[0,0,1,0,0],-1);
    P=dec2mat(lmis,xfeas,P)
    Q1=dec2mat(lmis,xfeas,Q1)
    Q2=dec2mat(lmis,xfeas,Q2)
    
    hH1=(hE+hL*hC)*inv(P)*Q1;
    hH2=(hE+hL*hC)*inv(P)*Q2;
    M=(hC*inv(P))';

    FFF=h1(n)*lmd1*yjy*jzjhx0';
    GGG=h2(n)*lmd2*yjy*jzjhx0';
    C11(n)=FFF(1,1);C12(n)=FFF(1,2);C13(n)=FFF(1,3);C14(n)=FFF(1,4);C15(n)=FFF(1,5);C16(n)=FFF(1,6);
    C21(n)=FFF(2,1);C22(n)=FFF(2,2);C23(n)=FFF(2,3);C24(n)=FFF(2,4);C25(n)=FFF(2,5);C26(n)=FFF(2,6);
    C31(n)=FFF(3,1);C32(n)=FFF(3,2);C33(n)=FFF(3,3);C34(n)=FFF(3,4);C35(n)=FFF(3,5);C36(n)=FFF(3,6);
    D11(n)=GGG(1,1);D12(n)=GGG(1,2);D13(n)=GGG(1,3);D14(n)=GGG(1,4);D15(n)=GGG(1,5);D16(n)=GGG(1,6);
    D21(n)=GGG(2,1);D22(n)=GGG(2,2);D23(n)=GGG(2,3);D24(n)=GGG(2,4);D25(n)=GGG(2,5);D26(n)=GGG(2,6);
    D31(n)=GGG(3,1);D32(n)=GGG(3,2);D33(n)=GGG(3,3);D34(n)=GGG(3,4);D35(n)=GGG(3,5);D36(n)=GGG(3,6);
    fff11=G.*C11/gamma(alpha);fff12=G.*C12/gamma(alpha);fff13=G.*C13/gamma(alpha);
    fff14=G.*C14/gamma(alpha);fff15=G.*C15/gamma(alpha);fff16=G.*C16/gamma(alpha);
    fff21=G.*C21/gamma(alpha);fff22=G.*C22/gamma(alpha);fff23=G.*C23/gamma(alpha);
    fff24=G.*C24/gamma(alpha);fff25=G.*C25/gamma(alpha);fff26=G.*C26/gamma(alpha);
    fff31=G.*C31/gamma(alpha);fff32=G.*C32/gamma(alpha);fff33=G.*C33/gamma(alpha);
    fff34=G.*C34/gamma(alpha);fff35=G.*C35/gamma(alpha);fff36=G.*C36/gamma(alpha);
    ggg11=G.*D11/gamma(alpha);ggg12=G.*D12/gamma(alpha);ggg13=G.*D13/gamma(alpha);
    ggg14=G.*D14/gamma(alpha);ggg15=G.*D15/gamma(alpha);ggg16=G.*D16/gamma(alpha);
    ggg21=G.*D21/gamma(alpha);ggg22=G.*D22/gamma(alpha);ggg23=G.*D23/gamma(alpha);
    ggg24=G.*D24/gamma(alpha);ggg25=G.*D25/gamma(alpha);ggg26=G.*D26/gamma(alpha);
    ggg31=G.*D31/gamma(alpha);ggg32=G.*D32/gamma(alpha);ggg33=G.*D33/gamma(alpha);
    ggg34=G.*D34/gamma(alpha);ggg35=G.*D35/gamma(alpha);ggg36=G.*D36/gamma(alpha);
    a11(n+1)=h*sum(fff11);a12(n+1)=h*sum(fff12);a13(n+1)=h*sum(fff13);
    a14(n+1)=h*sum(fff14);a15(n+1)=h*sum(fff15);a16(n+1)=h*sum(fff16);
    a21(n+1)=h*sum(fff21);a22(n+1)=h*sum(fff22);a23(n+1)=h*sum(fff23);
    a24(n+1)=h*sum(fff24);a25(n+1)=h*sum(fff25);a26(n+1)=h*sum(fff26);
    a31(n+1)=h*sum(fff31);a32(n+1)=h*sum(fff32);a33(n+1)=h*sum(fff33);
    a34(n+1)=h*sum(fff34);a35(n+1)=h*sum(fff35);a36(n+1)=h*sum(fff36);
    b11(n+1)=h*sum(ggg11);b12(n+1)=h*sum(ggg12);b13(n+1)=h*sum(ggg13);
    b14(n+1)=h*sum(ggg14);b15(n+1)=h*sum(ggg15);b16(n+1)=h*sum(ggg16);
    b21(n+1)=h*sum(ggg21);b22(n+1)=h*sum(ggg22);b23(n+1)=h*sum(ggg23);
    b24(n+1)=h*sum(ggg24);b25(n+1)=h*sum(ggg25);b26(n+1)=h*sum(ggg26);
    b31(n+1)=h*sum(ggg31);b32(n+1)=h*sum(ggg32);b33(n+1)=h*sum(ggg33);
    b34(n+1)=h*sum(ggg34);b35(n+1)=h*sum(ggg35);b36(n+1)=h*sum(ggg36);

    
 jhth1=[0+h*sum(fff11),0+h*sum(fff12),0+h*sum(fff13),0+h*sum(fff14),0+h*sum(fff15),0+h*sum(fff16);0+h*sum(fff21),0+h*sum(fff22),0+h*sum(fff23),0+h*sum(fff24),0+h*sum(fff25),0+h*sum(fff26);0+h*sum(fff31),0+h*sum(fff32),0+h*sum(fff33),0+h*sum(fff34),0+h*sum(fff35),0+h*sum(fff36)];
 jhth2=[0+h*sum(ggg11),0+h*sum(ggg12),0+h*sum(ggg13),0+h*sum(ggg14),0+h*sum(ggg15),0+h*sum(ggg16);0+h*sum(ggg21),0+h*sum(ggg22),0+h*sum(ggg23),0+h*sum(ggg24),0+h*sum(ggg25),0+h*sum(ggg26);0+h*sum(ggg31),0+h*sum(ggg32),0+h*sum(ggg33),0+h*sum(ggg34),0+h*sum(ggg35),0+h*sum(ggg36)];

    jhdA11=M*jhth1;
    jhdA12=M*jhth2;
    jhdA1=(hE+hL*hC)*jhdA11;
    jhdA2=(hE+hL*hC)*jhdA12;


    if t(n)>=2
    df1(n)=5*pi*sin(5*pi*(t(n)-5))*cos(5*pi*(t(n)-5));
    dff1=(t(n+1)-t(1:n)).^(-alpha).*df1/gamma(1-alpha);
    ddff1=h*sum(dff1);
    else
        ddff1=0;
    end
    if t(n)>=4
    df2(n)=-0.8*pi*sin(2*pi*(t(n)-4))*cos(2*pi*(t(n)-4));
    dff2=(t(n+1)-t(1:n)).^(-alpha).*df2/gamma(1-alpha);
    ddff2=h*sum(dff2);

    else
        ddff2=0;
    end
    df3(n)=1.8*pi*sin(3*pi*t(n))*cos(3*pi*t(n))-0.6*pi*sin(2*pi*t(n));
    dff3=(t(n+1)-t(1:n)).^(-alpha).*df3/gamma(1-alpha);
    ddff3=h*sum(dff3);
    HHH=[ddff1;ddff2;ddff3];
    CCC=h1(n)*C1*AAA+h2(n)*C2*AAA+HHH;
    DDD=hh*h1(n)*(hA01*jzjhx0+hB1*u(n)+hH1*yjy+jhdA1*jzjhx0+hN*[y1(n);y2(n);y3(n)]+hL*CCC)+hh*h2(n)*(hA02*jzjhx0+hB2*u(n)+hH2*yjy+jhdA2*jzjhx0+hN*[y1(n);y2(n);y3(n)]+hL*CCC);
    g1(n)=DDD(1);
    g2(n)=DDD(2);
    g3(n)=DDD(3);
    g4(n)=DDD(4);
    g5(n)=DDD(5);
    g6(n)=DDD(6);
    gg1=G.*g1/gamma(alpha);
    gg2=G.*g2/gamma(alpha);
    gg3=G.*g3/gamma(alpha);
    gg4=G.*g4/gamma(alpha);
    gg5=G.*g5/gamma(alpha);
    gg6=G.*g6/gamma(alpha);
    jx1(n+1)=jx1(1)+h*sum(gg1);
    jx2(n+1)=jx2(1)+h*sum(gg2);
    jx3(n+1)=jx3(1)+h*sum(gg3);
    jw01(n+1)=jw01(1)+h*sum(gg4);
    jw02(n+1)=jw02(1)+h*sum(gg5);
    jw03(n+1)=jw03(1)+h*sum(gg6);

     jzjw=[jw01(n);jw02(n);jw03(n)]-h1(n)*(C1-C)*[jx1(n);jx2(n);jx3(n)]-h2(n)*(C2-C)*[jx1(n);jx2(n);jx3(n)];
     jw1(n)=jzjw(1);
     jw2(n)=jzjw(2);
     jw3(n)=jzjw(3);
     e1(n)=x1(n)-jx1(n);
     e2(n)=x2(n)-jx2(n);
     e3(n)=x3(n)-jx3(n);
     e4(n)=w1(n)-jw1(n);
     e5(n)=w2(n)-jw2(n);
     e6(n)=w3(n)-jw3(n);


end

   
%plot(t,x1,'-',t,x2,'--',t,x3,':','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('State variables ')
%legend('x_1(t)','x_2(t)','x_3(t)')
%legend('boxoff')

%plot(t(1:n),y1,'-',t(1:n),y2,'--',t(1:n),y3,':','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Output variables ')
%legend('y_1(t)','y_2(t)','y_3(t)')
%legend('boxoff')


%plot(t,x1,'-',t,jx1,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $x_1$ and $\hat{x}_1$','interpreter','latex')
%legend('$x_1$','$\hat{x}_1$','interpreter','latex')
%legend('boxoff')

%plot(t,x2,'-',t,jx2,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $x_2$ and $\hat{x}_2$','interpreter','latex')
%legend('$x_2$','$\hat{x}_2$','interpreter','latex')
%legend('boxoff')

%plot(t,x3,'-',t,jx3,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $x_3$ and $\hat{x}_3$','interpreter','latex')
%legend('$x_3$','$\hat{x}_3$','interpreter','latex')
%legend('boxoff')

%plot(t(1:n),w1,'-',t(1:n),jw1,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $d_1$ and $\hat{d}_1$','interpreter','latex')
%legend('$d_1$','$\hat{d}_1$','interpreter','latex')
%legend('boxoff')


%plot(t(1:n),w2,'-',t(1:n),jw2,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $d_2$ and $\hat{d}_2$','interpreter','latex')
%legend('$d_2$','$\hat{d}_2$','interpreter','latex')
%legend('boxoff')

%plot(t(1:n),w3,'-',t(1:n),jw3,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $d_3$ and $\hat{d}_3$','interpreter','latex')
%legend('$d_3$','$\hat{d}_3$','interpreter','latex')
%legend('boxoff')
    
 %plot(t(1:n),e1,'-',t(1:n),e2,'-',t(1:n),e3,'-','LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_1','e_2','e_3')
 %legend('boxoff')

 %plot(t(1:n),e4,'-',t(1:n),e5,'-',t(1:n),e6,'-','LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_4','e_5','e_6')
 %legend('boxoff')


plot(t,x1,'-',t,jx1,'--','LineWidth',1.8)
xlabel('Time(second)')
ylabel('$x_1$ and $\hat{x}_1$','interpreter','latex')
legend('$x_1$','$\hat{x}_1$','interpreter','latex')
legend('boxoff')



% subplot(2,2,1)
% plot(t,x1,'-',t,jx1,'--','LineWidth',1.8)
% title('(a)')
% xlabel('Time(second)')
% ylabel('$x_1$ and $\hat{x}_1$','interpreter','latex')
% legend('$x_1$','$\hat{x}_1$','interpreter','latex')
% legend('boxoff')
% subplot(2,2,2)
% plot(t,x2,'-',t,jx2,'--','LineWidth',1.8)
% title('(b)')
% xlabel('Time(second)')
% ylabel('$x_2$ and $\hat{x}_2$','interpreter','latex')
% legend('$x_2$','$\hat{x}_2$','interpreter','latex')
% legend('boxoff')
% subplot(2,2,3)
% plot(t,x3,'-',t,jx3,'--','LineWidth',1.8)
% title('(c)')
% xlabel('Time(second)')
% ylabel('$x_3$ and $\hat{x}_3$','interpreter','latex')
% legend('$x_3$','$\hat{x}_3$','interpreter','latex')
% legend('boxoff')
% subplot(2,2,4)
% plot(t(1:n),e1,'-',t(1:n),e2,'-',t(1:n),e3,'-','LineWidth',1.8)
% title('(d)')
% %ylim([-1 1])
% xlabel('Time(second)')
% ylabel('Tracking errors')
% legend('e_1','e_2','e_3')
% legend('boxoff')

% plot(t(1:n),e1,'-',t(1:n),e2,'-',t(1:n),e3,'-','LineWidth',1.8)
% ylim([-2 2])
% xlabel('Time(second)')
% ylabel('Tracking errors')
% legend('e_1','e_2','e_3')
% legend('boxoff')






%subplot(2,2,1)
%plot(t(1:n),w1,'-',t(1:n),jw1,'--','LineWidth',1.8)
%title('(a)')
%xlabel('Time(second)')
%ylabel('$d_1$ and $\hat{d}_1$','interpreter','latex')
%legend('$d_1$','$\hat{d}_1$','interpreter','latex')
%legend('boxoff')
%subplot(2,2,2)
%plot(t(1:n),w2,'-',t(1:n),jw2,'--','LineWidth',1.8)
%title('(b)')
%xlabel('Time(second)')
%ylabel('$d_2$ and $\hat{d}_2$','interpreter','latex')
%legend('$d_2$','$\hat{d}_2$','interpreter','latex')
%legend('boxoff')
%subplot(2,2,3)
%plot(t(1:n),w3,'-',t(1:n),jw3,'--','LineWidth',1.8)
%title('(c)')
%xlabel('Time(second)')
%ylabel('$d_3$ and $\hat{d}_3$','interpreter','latex')
%legend('$d_3$','$\hat{d}_3$','interpreter','latex')
%legend('boxoff')
%subplot(2,2,4)
%plot(t(1:n),e4,'-',t(1:n),e5,'-',t(1:n),e6,'-','LineWidth',1.8)
%title('(d)')
%ylim([-1 1])
%xlabel('Time(second)')
%ylabel('Tracking errors')
%legend('e_4','e_5','e_6')
%legend('boxoff')

















