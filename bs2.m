clear all
clc 
x1=0.35,x2=0.55,x3=-0.25,x4=0.65;
jx1=0,jx2=0,jx3=0,jx4=0;
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
    y(n)=BBB;

    fg1(n)=jx2(n)+8*(x1(n)+jzw-jx1(n));
    abc1=G.*fg1/gamma(alpha);
    jx1(n+1)=jx1(1)+h*sum(abc1);
    fg2(n)=jx3(n)-10*sin(jx1(n))-10*jx1(n)+24*(x1(n)+jzw-jx1(n));
    abc2=G.*fg2/gamma(alpha);
    jx2(n+1)=jx2(1)+h*sum(abc2);
    fg3(n)=jx4(n)+24*(x1(n)+jzw+jx1(n));
    abc3=G.*fg3/gamma(alpha);
    jx3(n+1)=jx3(1)+h*sum(abc3);
    fg4(n)=10*jx1(n)-2*jx3(n)+4*u(n)+16*(x1(n)+jzw-jx1(n));
    abc4=G.*fg4/gamma(alpha);
    jx4(n+1)=jx4(1)+h*sum(abc4);





     e1(n)=x1(n)-jx1(n);
     e2(n)=x2(n)-jx2(n);
     e3(n)=x3(n)-jx3(n);
     e4(n)=x4(n)-jx4(n);

     


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

%plot(t(1:n),w3,'-',t(1:n),jw3,'--','LineWidth',1.8)
%xlabel('Time(second)')
%ylabel('Tracking tractories of $d_3$ and $\hat{d}_3$','interpreter','latex')
%legend('$d_3$','$\hat{d}_3$','interpreter','latex')
%legend('boxoff')
    
 %plot(t(1:n),e1,'-',t(1:n),e2,'-',t(1:n),e3,'-',t(1:n),e4,'LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_1','e_2','e_3','e_4')
 %legend('boxoff')

 %plot(t(1:n),e4,'-',t(1:n),e5,'-',t(1:n),e6,'-','LineWidth',1.8)
 %ylim([-1 1])
 %xlabel('Time(second)')
 %ylabel('Tracking errors')
 %legend('e_4','e_5','e_6')
 %legend('boxoff')
































