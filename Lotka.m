function [xsol ysol] = Lotka(xn,yn,alpha,beta,gamma,delta,n,ntirage,npasX,npasY,Cbio,sigmaB,sigmaH)
h=1/n;T=n*10;% Model parameters


k1x=xn*(alpha-beta*yn);
k1y=yn*(delta*xn-gamma);
k2x=(xn+h/2*k1x)*(alpha-beta*(yn+h/2*k1y));
k2y=(yn+h/2*k1y)*(delta*(xn+h/2*k1x)-gamma);
k3x=(xn+h/2*k2x)*(alpha-beta*(yn+h/2*k2y));
k3y=(yn+h/2*k2y)*(delta*(xn+h/2*k2x)-gamma);
k4x=(xn+h*k3x)*(alpha-beta*(yn+h*k3y));
k4y=(yn+h*k3y)*(delta*(xn+h*k3x)-gamma);
xn=xn+h/6*(k1x+2*k2x+2*k3x+k4x);
yn=yn+h/6*(k1y+2*k2y+2*k3y+k4y);

xn=xn+randn*sigmaB;
yn=yn+randn*sigmaH;
if xn<0
    xn=0;
end

if yn<0
    yn=0;
end      
xsol=xn;ysol=yn;