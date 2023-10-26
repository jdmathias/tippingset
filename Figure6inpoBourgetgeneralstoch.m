close all;clear all;
alpha=0.0495;beta=-23.693;s=2.1476;h=0.12;r=367.04;q=2.222;m=96.85;
T=12;dt=1;
ntirage=1000;npasX=200;npasY=200;%1000 x 3
xm=180;ym=2600;pasX=xm/npasX;pasY=ym/npasY;[X,Y] = meshgrid(pasX/2:pasX:(xm-pasX/2), pasY/2:pasY:(ym-pasY/2));
Cbio=42;
%%%%%%%%% 
xdata=[118.8  111.6 86.4 79.2 72 72 61.2 57.6 50.4 39.6 39.6 36 28.8]
ydata=[1027 888 1004 1238 1214 929 1031 1079 1385 1467 1231 1332 1181]
mur=mean(ydata);
sigmar=std(ydata);%% attention je trouve 176,84
betat=0.75;gammat=0.1;% MFET threshold; 
limu=2500;
%%%%%%%%% loop for each point of the grid
for ix=1:npasX % loop on biomass values
    ix
    for jy=1:npasY % loop on effort values
        for k=1:ntirage        
            x(1)=X(jy,ix);y(1)=Y(jy,ix);
            critere=0;critereu=0;cini(ix,jy)=0;inD(ix,jy,k)=0;inU(ix,jy,k)=0;
            if x(1)>Cbio 
                 critereu=1;
                 inU(ix,jy,k)=1;
                 cini(ix,jy)=-1;
            else        
                 critere=1; 
                 inD(ix,jy,k)=1;  
            end        
            for t=2:T %%%%%% time loop
                Xp=x(t-1);Yp=y(t-1);
                for tat=1:10
                     Xp=Xp+(alpha*Yp+beta-(s+h)*Xp+r*Xp^q/(Xp^q+m^q))*dt/10;%small time step for taking into account non linearities of the dynamics
                end
                x(t)=Xp;
                if x(t)<0
                    x(t)=0;
                end
                y(t)=randn*sigmar+mur;
                if y(t)<0
                    y(t)=0;
                end

                if x(t)>Cbio && critereu==0
                    critereu=1;
                    FET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if x(t)<Cbio && critere==0
                    critere=1;
                    FETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end 
                if x(t)>Cbio 
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
            
            
            for t=2:round(gammat*T*1.1+2) %%%%%% time loop
                Xp=x(t-1);Yp=y(t-1);
                for tat=1:10
                     Xp=Xp+(alpha*Yp+beta-(s+h)*Xp+r*Xp^q/(Xp^q+m^q))*dt/10;%small time step for taking into account non linearities of the dynamics
                end
                x(t)=Xp;
                if x(t)<pasX/2
                    x(t)=pasX/2;
                end
                y(t)=randn*sigmar+mur;
                if y(t)<pasY/2
                    y(t)=pasY/2;
                end

                %%%%%%%%% MFT
                if inpolygon(y(t),x(t),[xMFETd1 xMFETd2(end:-1:1)],[VMFETdmin VMFETdmax(end:-1:1)])==1 && critereMFET==1
                    critereMFET=0;
                    FATMFET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0   
                end
                if  inpolygon(y(t),x(t),[xMFETu1 xMFETu2(end:-1:1)],[VMFETumin VMFETumax(end:-1:1)])==1 && critereuMFET==1
                    critereuMFET=0;
                    FATMFETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if  x(t)>Cbio && y(t) >limu && critereuMFET==1
                    critereuMFET=0;
                    FATMFETu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                
                if t==round(gammat*T*1.1+2) && critereMFET==1
                    critereMFET=0;
                    FATMFET(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==round(gammat*T*1.1+2) && critereuMFET==1
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
                if  x(t)>Cbio && y(t) >limu && critereuMPT==1
                    critereuMPT=0;
                    FATMPTu(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==round(gammat*T*1.1+2) && critereMPT==1
                    critereMPT=0;
                    FATMPT(ix,jy,k)=t-1;%we remove 1 because we start at time t=1 instead of t=0 
                end
                if t==round(gammat*T*1.1+2) && critereuMPT==1
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

xeq=1:1:180;
yeq=-(+beta-(s+h)*xeq+r*xeq.^q./(xeq.^q+m^q))/alpha;

figure; 
M=MPT*0+NaN;
MPTtest=MPT;MPTtest(MPTtest<=betat)=NaN;
MPTutest=MPTu;MPTutest(MPTutest<=betat)=NaN;
MFATMPTtest=MFATMPT;MFATMPTtest(MFATMPTtest>gammat)=NaN;
MFATMPTutest=MFATMPTu;MFATMPTutest(MFATMPTtest>gammat)=NaN;
M(MPTtest>betat)=1;M(MPTutest>betat)=0;
M(MFATMPTtest<=gammat)=-0.5;M(MFATMPTutest<=gammat)=0.5;
clims = [-1 1];imagesc(Y(:,1),X(1,:),M,clims);set(gca,'YDir','normal')
colormap([1 1 1;0.3010 0.7450 0.9330 ; 1 0 0; 1 0.8 0.8; 0 0 1] )
hold on
plot(yeq,xeq,'k','Linewidth',3)
plot(ydata,xdata,'g','Linewidth',3)
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
hold on
plot(yeq,xeq,'k','Linewidth',3)
plot(ydata,xdata,'g','Linewidth',3)
ylim([0 xm]);xlim([0 ym])