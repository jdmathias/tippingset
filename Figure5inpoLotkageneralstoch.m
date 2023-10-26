close all;clear all;
alpha=1;beta=0.75;gamma=1;delta=0.75;n=100;h=1/n;T=10*n;dt=1;% Model parameters:n=100,ntirage=1000,npasX=npasy=500 en det et 200 en stoch
ntirage=1;npasX=500;npasY=500;%nbre tirage; pas de grille. Pour les simulations: 1000 tirages et 200 pasx et 200 pasy
xm=8;ym=9;pasX=xm/npasX;pasY=ym/npasY;[X,Y] = meshgrid(0+pasX/2:pasX:(xm-pasX/2), pasY/2:pasY:(ym-pasY/2));
Cbio=1;sigmaB=0.005*0;sigmaH=0.005*0;
betat=0.6;gammat=0.1;%
limu=8.8;
%%%%%%%%% loop for each point of the grid
for ix=1:npasX % loop on biomass values
    ix
    for jy=1:npasY % loop on effort values
        for k=1:ntirage        
            x(1)=X(jy,ix);y(1)=Y(jy,ix);
            critere=0;critereu=0;cini(ix,jy)=0;inD(ix,jy,k)=0;inU(ix,jy,k)=0;
            if x(1)<Cbio 
                 critereu=1;
                 inU(ix,jy,k)=1;
                 cini(ix,jy)=-1;
            else        
                 critere=1; 
                 inD(ix,jy,k)=1;  
            end        
            for t=2:T %%%%%% time loop
                yn=y(t-1);xn=x(t-1);
                [x(t) y(t)]=Lotka(xn,yn,alpha,beta,gamma,delta,n,ntirage,npasX,npasY,Cbio,sigmaB,sigmaH);
                if x(t)<Cbio && critereu==0
                    critereu=1;
                    FET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if x(t)>Cbio && critere==0
                    critere=1;
                    FETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end 
                if x(t)<Cbio 
                    inU(ix,jy,k)=inU(ix,jy,k)+1;
                else
                    inD(ix,jy,k)=inD(ix,jy,k)+1;
                end 
                if t==T && critereu==0
                    critereu=1;
                    FET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==T && critere==0
                    critere=1;
                    FETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
            end
        end
        if cini(ix,jy)==0
            MFET(ix,jy)=median(FET(ix,jy,:))/T;
            MFETu(ix,jy)=NaN;
            MPT(ix,jy)=mean(inD(ix,jy,:))/T;
            MPTu(ix,jy)=0;
        else
            MFETu(ix,jy)=median(FETu(ix,jy,:))/T;
            MFET(ix,jy)=NaN;
            MPT(ix,jy)=0;
            MPTu(ix,jy)=mean(inU(ix,jy,:))/T;
        end    
    end
end
DMFET=X'*0+1*(MFET>betat);%ensemble desirable
UMFET=X'*0+1*(MFETu>betat);
DMPT=X'*0+1*(MPT>betat);
UMPT=X'*0+1*(MPTu>betat);

%%%%%%%% Set DMFET
k=0;xMFETd1=0;xMFETd2=0;VMFETdmin=0;VMFETdmax=0;
for j=1:npasY
    crmin=0;
    for i=1:npasX
        if crmin==0 && DMFET(i,j)==1
            k=k+1;
            VMFETdmin(k)=X(1,i);%i;
            xMFETd1(k)=Y(j,1);%j
            crmin=1;
        end
        if DMFET(i,j)==1
            VMFETdmax(k)=X(1,i);
            xMFETd2(k)=Y(j,1);
        end
    end
end    
%%%%%%%% Set UMFET
k=0;xMFETu1=0;xMFETu2=0;VMFETumin=0;VMFETumax=0;
for j=1:npasY
    crmin=0;
    for i=1:npasX
        if crmin==0 && UMFET(i,j)==1
            k=k+1;
            VMFETumin(k)=X(1,i);%i;
            xMFETu1(k)=Y(j,1);%j
            crmin=1;
        end
        if UMFET(i,j)==1
            VMFETumax(k)=X(1,i);
            xMFETu2(k)=Y(j,1);
        end
    end
end    

%%%%%%%% Set DMPT
k=0;xMPTd1=0;xMPTd2=0;VMPTdmin=0;VMPTdmax=0;
for j=1:npasY
    crmin=0;
    for i=1:npasX
        if crmin==0 && DMPT(i,j)==1
            k=k+1;
            VMPTdmin(k)=X(1,i);%i;
            xMPTd1(k)=Y(j,1);%j
            crmin=1;
        end
        if DMPT(i,j)==1
            VMPTdmax(k)=X(1,i);
            xMPTd2(k)=Y(j,1);
        end
    end
end    
%%%%%%%% Set UMPT
k=0;xMPTu1=0;xMPTu2=0;VMPTumin=0;VMPTumax=0;
for j=1:npasY
    crmin=0;
    for i=1:npasX
        if crmin==0 && UMPT(i,j)==1
            k=k+1;
            VMPTumin(k)=X(1,i);%i;
            xMPTu1(k)=Y(j,1);%j
            crmin=1;
        end
        if UMPT(i,j)==1
            VMPTumax(k)=X(1,i);
            xMPTu2(k)=Y(j,1);
        end
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reachability
%%%%%%%%% loop for each point of the grid
for ix=1:npasX % loop on biomass values
    ix
    for jy=1:npasY % loop on effort values
        for k=1:ntirage
            x(1)=X(jy,ix);y(1)=Y(jy,ix);
            critereMFET=0;critereuMFET=0;
            critereMPT=0;critereuMPT=0;
            ciniMFET(ix,jy)=0;ciniuMFET(ix,jy)=0;
            ciniMPT(ix,jy)=0;ciniuMPT(ix,jy)=0;
            if DMFET(ix,jy)==0 
                 critereMFET=1;
                 ciniMFET(ix,jy)=-1;
            end
            if UMFET(ix,jy)==0 
                 critereuMFET=1;
                 ciniuMFET(ix,jy)=-1;
            end
            if DMPT(ix,jy)==0 
                 critereMPT=1;
                 ciniMPT(ix,jy)=-1;
            end
            if UMPT(ix,jy)==0 
                 critereuMPT=1;
                 ciniuMPT(ix,jy)=-1;
            end
            
            for t=2:round(gammat*T*1.1) %%%%%% time loop
                yn=y(t-1);xn=x(t-1);
                [x(t) y(t)]=Lotka(xn,yn,alpha,beta,gamma,delta,n,ntirage,npasX,npasY,Cbio,sigmaB,sigmaH);
                %%%%%%%%% MFT
                if inpolygon(y(t),x(t),[xMFETd1 xMFETd2(end:-1:1)],[VMFETdmin VMFETdmax(end:-1:1)])==1 && critereMFET==1
                    critereMFET=0;
                    FATMFET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0   
                end
                if  inpolygon(y(t),x(t),[xMFETu1 xMFETu2(end:-1:1)],[VMFETumin VMFETumax(end:-1:1)])==1 && critereuMFET==1
                    critereuMFET=0;
                    FATMFETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if  x(t)<Cbio && y(t) >limu && critereuMFET==1
                    critereuMFET=0;
                    FATMFETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                
                if t==round(gammat*T*1.1) && critereMFET==1
                    critereMFET=0;
                    FATMFET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==round(gammat*T*1.1) && critereuMFET==1
                    critereuMFET=0;
                    FATMFETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                %%%%%%%%%%%%%% MPT
                if inpolygon(y(t),x(t),[xMPTd1 xMPTd2(end:-1:1)],[VMPTdmin VMPTdmax(end:-1:1)])==1 && critereMPT==1
                    critereMPT=0;
                    FATMPT(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0   
                end
                if  inpolygon(y(t),x(t),[xMPTu1 xMPTu2(end:-1:1)],[VMPTumin VMPTumax(end:-1:1)])==1 && critereuMPT==1
                    critereuMPT=0;
                    FATMPTu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if  x(t)<Cbio && y(t) >limu && critereuMPT==1
                    critereuMPT=0;
                    FATMPTu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end                
                if t==round(gammat*T*1.1) && critereMPT==1
                    critereMPT=0;
                    FATMPT(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==round(gammat*T*1.1) && critereuMPT==1
                    critereuMPT=0;
                    FATMPTu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end  
            end   
        end
        
        if ciniMFET(ix,jy)==-1
            MFATMFET(ix,jy)=median(FATMFET(ix,jy,:))/T;
        else
            MFATMFET(ix,jy)=NaN;
        end
        if ciniuMFET(ix,jy)==-1
            MFATMFETu(ix,jy)=median(FATMFETu(ix,jy,:))/T;
        else
            MFATMFETu(ix,jy)=NaN;
        end
        
         if ciniMPT(ix,jy)==-1
            MFATMPT(ix,jy)=median(FATMPT(ix,jy,:))/T;
        else
            MFATMPT(ix,jy)=NaN;
        end
        if ciniuMPT(ix,jy)==-1
            MFATMPTu(ix,jy)=median(FATMPTu(ix,jy,:))/T;
        else
            MFATMPTu(ix,jy)=NaN;
        end
        
    end
end

RDMFET=X'*0+1*(MFATMFET<=gammat);
RUMFET=X'*0+1*(MFATMFETu<=gammat);

RDMPT=X'*0+1*(MFATMPT<=gammat);
RUMPT=X'*0+1*(MFATMPTu<=gammat);

figure; 
M=MPT*0+NaN;
MPTtest=MPT;MPTtest(MPTtest<=betat)=NaN;%MPTtest(MPTtest>betat)=1;
MPTutest=MPTu;MPTutest(MPTutest<=betat)=NaN;%MPTutest(MPTutest>betat)=1;
MFATMPTtest=MFATMPT;MFATMPTtest(MFATMPTtest>gammat)=NaN;%MPTtest(MPTtest>betat)=1;
MFATMPTutest=MFATMPTu;MFATMPTutest(MFATMPTtest>gammat)=NaN;
M(MPTtest>betat)=1;M(MPTutest>betat)=0;
M(MFATMPTtest<=gammat)=-0.5;M(MFATMPTutest<=gammat)=0.5;
clims = [-1 1];imagesc(Y(:,1),X(1,:),M,clims);set(gca,'YDir','normal')
colormap([1 1 1;0.3010 0.7450 0.9330 ; 1 0 0; 1 0.8 0.8; 0 0 1])
ylim([0 xm]);xlim([0 ym])

figure;
M=MFET*0+NaN;
MFETtest=MFET;MFETtest(MFETtest<=betat)=NaN;
MFETutest=MFETu;MFETutest(MFETutest<=betat)=NaN;
MFATMFETtest=MFATMFET;MFATMFETtest(MFATMFETtest>gammat)=NaN;
MFATMFETutest=MFATMFETu;MFATMFETutest(MFATMFETtest>gammat)=NaN;
M(MFETtest>betat)=1;M(MFETutest>betat)=0;
M(MFATMFETtest<=gammat)=-0.5;M(MFATMFETutest<=gammat)=0.5;
clims = [-1 1];imagesc(Y(:,1),X(1,:),M,clims);set(gca,'YDir','normal')
colormap([1 1 1;0.3010 0.7450 0.9330 ; 1 0 0; 1 0.8 0.8; 0 0 1] )
ylim([0 xm]);xlim([0 ym])

x(1)=1.5;y(1)=2.5;
for t=2:T %%%%%% time loop
        yn=y(t-1);xn=x(t-1);
        [x(t) y(t)]=Lotka(xn,yn,alpha,beta,gamma,delta,n,ntirage,npasX,npasY,Cbio,sigmaB,sigmaH);
end
hold on; plot(y,x);


