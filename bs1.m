clear all
clc 
x1=0.35,x2=0.55,x3=-0.25,x4=0.65;
jx1=0,jx2=0,jx3=0,jx4=0;
jw0=0; 
jy=0;  
t=0;
h=0.0001;
alpha=0.95;
A1=[0,1,0,0;-20,0,1,0;0,0,0,1;10,0,-2,0];
A2=[0,1,0,0;-10,0,1,0;0,0,0,1;10,0,-2,0];
B1=[0;0;0;4];
B2=[0;0;0;4];
C1=[1,0,0,0];%C2不等C1;
C2=[-1,0,0,0];
C=C1;
lmd1=1;
lmd2=2;
dlA1=[0.0708,0,0,0;0,0.05,0,0;0,0,0.03,0;0,0,0,0.01];
dlA2=[0.2326,0,0,0;0,0.2,0,0;0,0,0.03,0;0,0,0,0.01];
epsl1=norm(dlA1);
epsl2=norm(dlA2);
gam1=0.01;
gam2=0.02;
hE=[1,0,0,0,0;0,1,0,0,0;0,0,1,0,0;0,0,0,1,0;0,0,0,0,0];
hA1=[0,1,0,0,0;-20,0,1,0,0;0,0,0,1,0;10,0,-2,0,0;0,0,0,0,-1];
hA2=[0,1,0,0,0;-10,0,1,0,0;0,0,0,1,0;10,0,-2,0,0;0,0,0,0,-1];
hB1=[0;0;0;4;0];hB2=[0;0;0;4;0];
hN=[0;0;0;0;1];
hC=[1,0,0,0,1];
C0=[1,0,0,0,0];
hL=[0;0;0;0;1];
hh=inv(hE+hL*hC);
epsl11=norm(hh)*epsl1;
epsl12=norm(hh)*epsl2;
p1=epsl11/(norm(hC)*gam1);
p2=epsl12/(norm(hC)*gam2);
a11=0,a12=0,a13=0,a14=0,a15=0,a16=0,a21=0,a22=0,a23=0,a24=0,a25=0,a26=0,a31=0,a32=0,a33=0,a34=0,a35=0,a36=0;
b11=0,b12=0,b13=0,b14=0,b15=0,b16=0,b21=0,b22=0,b23=0,b24=0,b25=0,b26=0,b31=0,b32=0,b33=0,b34=0,b35=0,b36=0;
jhth1=[a11,a12,a13,a14,a15,a16;a21,a22,a23,a24,a25,a26;a31,a32,a33,a34,a35,a36];
jhth2=[b11,b12,b13,b14,b15,b16;b21,b22,b23,b24,b25,b26;b31,b32,b33,b34,b35,b36];


for n=1:199990
    t(n+1)=t(n)+h;
    u(n)=sin(t(n));
    h1(n)=0.5*(1+0.1*(4-x1(n)));
    h2(n)=0.5*(1-0.1*(4-x1(n)));
    dA1=dlA1*cos(t(n));
    dA2=dlA2*cos(t(n));
    hdA1=[0.0708,0,0,0,0;0,0.05,0,0,0;0,0,0.03,0,0;0,0,0,0.01,0;0,0,0,0,0]*cos(t(n));
    hdA2=[0.2326,0,0,0,0;0,0.2,0,0,0;0,0,0.03,0,0;0,0,0,0.01,0;0,0,0,0,0]*cos(t(n));
    jz=[x1(n);x2(n);x3(n);x4(n)];
    AAA=h1(n)*((A1+dA1)*jz+B1*u(n))+h2(n)*((A2+dA2)*jz+B2*u(n));
    f1(n)=AAA(1);
    f2(n)=AAA(2);
    f3(n)=AAA(3);
    f4(n)=AAA(4);
    G=(t(n+1)-t(1:n)).^(alpha-1);
    ff1=G.*f1/gamma(alpha);
    ff2=G.*f2/gamma(alpha);
    ff3=G.*f3/gamma(alpha);
    ff4=G.*f4/gamma(alpha);
    x1(n+1)=x1(1)+h*sum(ff1);
    x2(n+1)=x2(1)+h*sum(ff2);
    x3(n+1)=x3(1)+h*sum(ff3);
    x4(n+1)=x4(1)+h*sum(ff4);
    if t(n)>=3&&t(n)<=20
        w(n)=0.8*(sin(0.2*t(n)));
    else
        w(n)=0;
    end
    jzw=w(n);
    BBB=h1(n)*C1*jz+h2(n)*C2*jz+jzw;
    EEE=h1(n)*(C1-C)*jz+h2(n)*(C2-C)*jz+jzw;
    w0(n)=EEE;
    y(n)=BBB;
    yy=y(n);
    jyy=hC*[jx1(n);jx2(n);jx3(n);jx4(n);jw0(n)];
    hA01=hA1-hN*C0;
    hA02=hA2-hN*C0;
    jzjhx0=[jx1(n);jx2(n);jx3(n);jx4(n);jw0(n)];
    yjy=yy-jyy;
    hA11=hh*hA01;
    hA12=hh*hA02;

    setlmis([]);
    P=lmivar(1,[5,1]);
    Q1=lmivar(2,[5,1]);
    Q2=lmivar(2,[5,1]);
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
    C11(n)=FFF(1,1);C12(n)=FFF(1,2);C13(n)=FFF(1,3);C14(n)=FFF(1,4);C15(n)=FFF(1,5);
    D11(n)=GGG(1,1);D12(n)=GGG(1,2);D13(n)=GGG(1,3);D14(n)=GGG(1,4);D15(n)=GGG(1,5);
    fff11=G.*C11/gamma(alpha);fff12=G.*C12/gamma(alpha);fff13=G.*C13/gamma(alpha);
    fff14=G.*C14/gamma(alpha);fff15=G.*C15/gamma(alpha);
    ggg11=G.*D11/gamma(alpha);ggg12=G.*D12/gamma(alpha);ggg13=G.*D13/gamma(alpha);
    ggg14=G.*D14/gamma(alpha);ggg15=G.*D15/gamma(alpha);
    a11(n+1)=h*sum(fff11);a12(n+1)=h*sum(fff12);a13(n+1)=h*sum(fff13);
    a14(n+1)=h*sum(fff14);a15(n+1)=h*sum(fff15);
    b11(n+1)=h*sum(ggg11);b12(n+1)=h*sum(ggg12);b13(n+1)=h*sum(ggg13);
    b14(n+1)=h*sum(ggg14);b15(n+1)=h*sum(ggg15);

    
 jhth1=[0+h*sum(fff11),0+h*sum(fff12),0+h*sum(fff13),0+h*sum(fff14),0+h*sum(fff15)];
 jhth2=[0+h*sum(ggg11),0+h*sum(ggg12),0+h*sum(ggg13),0+h*sum(ggg14),0+h*sum(ggg15)];

    jhdA11=M*jhth1;
    jhdA12=M*jhth2;
    jhdA1=(hE+hL*hC)*jhdA11;
    jhdA2=(hE+hL*hC)*jhdA12;


    if t(n)>=3&&t(n)<=20
    df1(n)=0.16*cos(0.2*t(n));
    dff1=(t(n+1)-t(1:n)).^(-alpha).*df1/gamma(1-alpha);
    ddff1=h*sum(dff1);
    else
        ddff1=0;
    end
    HHH=ddff1;
    CCC=h1(n)*C1*AAA+h2(n)*C2*AAA+HHH;
    DDD=hh*h1(n)*(hA01*jzjhx0+hB1*u(n)+hH1*yjy+jhdA1*jzjhx0+hN*y(n)+hL*CCC)+hh*h2(n)*(hA02*jzjhx0+hB2*u(n)+hH2*yjy+jhdA2*jzjhx0+hN*y(n)+hL*CCC);
    g1(n)=DDD(1);
    g2(n)=DDD(2);
    g3(n)=DDD(3);
    g4(n)=DDD(4);
    g5(n)=DDD(5);
    gg1=G.*g1/gamma(alpha);
    gg2=G.*g2/gamma(alpha);
    gg3=G.*g3/gamma(alpha);
    gg4=G.*g4/gamma(alpha);
    gg5=G.*g5/gamma(alpha);
    jx1(n+1)=jx1(1)+h*sum(gg1);
    jx2(n+1)=jx2(1)+h*sum(gg2);
    jx3(n+1)=jx3(1)+h*sum(gg3);
    jx4(n+1)=jx4(1)+h*sum(gg4);
    jw0(n+1)=jw0(1)+h*sum(gg5);

     jzjw=jw0(n)-h1(n)*(C1-C)*[jx1(n);jx2(n);jx3(n);jx4(n)]-h2(n)*(C2-C)*[jx1(n);jx2(n);jx3(n);jx4(n)];
     jw(n)=jzjw;
     e1(n)=x1(n)-jx1(n);
     e2(n)=x2(n)-jx2(n);
     e3(n)=x3(n)-jx3(n);
     e4(n)=x4(n)-jx4(n);
     e5(n)=w(n)-jw(n);
     


end

   
%plot(t,x1,'-',t,x2,'--',t,x3,':','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('State variables ')
%legend('x_1(t)','x_2(t)','x_3(t)')
%legend('boxoff')

plot(t,x1,'-',t,jx1,'--','LineWidth',1.8)
xlabel('Time(second)')
ylabel('Tracking tractories of $x_1$ and $\hat{x}_1$','interpreter','latex')
legend('$x_1$','$\hat{x}_1$','interpreter','latex')
legend('boxoff')

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

%plot(t,x4,'-',t,jx4,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $x_4$ and $\hat{x}_4$','interpreter','latex')
%legend('$x_4$','$\hat{x}_4$','interpreter','latex')
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

%plot(t(1:n),w,'-',t(1:n),jw,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $d$ and $\hat{d}$','interpreter','latex')
%legend('$d$','$\hat{d}$','interpreter','latex')
%legend('boxoff')
    
 %plot(t(1:n),e1,'-',t(1:n),e2,'-',t(1:n),e3,'-','LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_1','e_2','e_3')
 %legend('boxoff')

 %plot(t(1:n),e4,'-',t(1:n),e5,'-','LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_4','e_5')
 %legend('boxoff')
































