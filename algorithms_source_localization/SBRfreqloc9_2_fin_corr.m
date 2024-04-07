function [I,orr,ShatI,GOF,XhatI] = SBRfreqloc9_2_fin_corr(X,K,lambda)
% SBRFREQLOC8 implements a block-SBR regression under rank 1 constranints
% V=data matrix, one signal per row
% L=lead-field
% expGOF=expected GOF
% indloc=source indices (positions)
% orr=source orientations
% Stilde=source amplitudes
% GOF=GOF

% SBR step sur Iold (avant addition)
% ligne 174: avec actualisation de lambda (CSBR)


normX=norm(X,'fro');

crit=Inf;
GOF=0;
ndip=size(K,2)/3;

I=[];
ShatI=[];
Jblk=[];

stop=0;
KI=[];
err_rec=Inf;

while ~stop
    
    % pour SBR sur old et pour savoir quoi sortir quand on converge
    ShatIold=ShatI;
    Jblkold=Jblk;
    KIold=KI;
    Iold=I;
    GOFold=GOF;
    
    % pour upper bound
    err=Inf(1,ndip);
    err2=Inf(1,ndip);
    I2=I;
    Jblk2=Jblk;
    KI2=KI;
    
    ndipcur=length(I)+1;
    
    for idip=setdiff([1:ndip],I)%
        ki=K(:,3*idip-2:3*idip);
        Shati=pinv([KI ki])*X;
        Jblki=zeros(3*ndipcur,ndipcur);
        for ic=1:ndipcur,
            [u,~,~]=svd(Shati(3*ic-2:3*ic,:));
            Jblki(3*ic-2:3*ic,ic)=u(:,1);
        end
        Xhati=[KI ki]*Jblki*pinv([KI ki]*Jblki)*X;
        err(idip)=norm(X-Xhati,'fro');
    end
    [errm,inddip1]=min(err);
    I=[I, inddip1];
    ki=K(:,3*inddip1-2:3*inddip1);
    KI=[KI ki];
    JI=pinv(KI)*X;
    Jblk=zeros(3*ndipcur,ndipcur);
    for ic=1:ndipcur,
        [U,~,~]=svd(JI(3*ic-2:3*ic,:));
        Jblk(3*ic-2:3*ic,ic)=U(:,1);
    end
    ShatI=pinv(KI*Jblk)*X;
    XhatI=KI*Jblk*ShatI; % these are the estimates after adding using normal path
    errm=norm(X-XhatI,'fro');
    
    %%
    
    for idip=setdiff([1:ndip],I2)%
        ki=K(:,3*idip-2:3*idip);
        Shat=pinv([KI2*Jblk2 ki])*X;
        [u,~,~]=svd(Shat(end-2:end,:));
        KI2t=[KI2 ki];
        Jblk2t=[Jblk2 zeros(3*(ndipcur-1),1); zeros(3,ndipcur-1) u(:,1)];
        Shati2t=pinv(KI2t*Jblk2t)*X;
        Xhati2=KI2t*Jblk2t*Shati2t;
        err2(idip)=norm(X-Xhati2,'fro');
    end
    [errm2,inddip2]=min(err2);
    I2=[I2 inddip2];
    ki=K(:,3*inddip2-2:3*inddip2);
    Shat=pinv([KI2*Jblk2 ki])*X;
    [u,~,~]=svd(Shat(end-2:end,:));
    KI2=[KI2 ki];
    Jblk2=[Jblk2 zeros(3*(ndipcur-1),1); zeros(3,ndipcur-1) u(:,1)];
    ShatI2=pinv(KI2*Jblk2)*X;
    XhatI2=KI2*Jblk2*ShatI2;
    errm2=norm(X-XhatI2,'fro');
    if errm2<errm, % these are the estimates after adding using alternate path (upper bounded)
        I=I2;
        Jblk=Jblk2;
        ShatI=ShatI2;
        KI=KI2;
        XhatI=XhatI2;
    end
    err_rec_old=err_rec;
    err_rec=norm(X-XhatI,'fro');
    GOF=1-(err_rec/normX)^2;
    
    critadd=err_rec/normX+lambda*length(I);
    %% SBR step
    critsub=Inf;
    if length(Iold)>1
        critsubv=Inf(1,length(Iold));
        for idip=1:length(Iold)
            indcur=Iold([1:idip-1 idip+1:length(Iold)]);
            iic=[1:idip-1 idip+1:length(Iold)];
            tmp=[3*iic-2;3*iic-1;3*iic];
            KIcur=KIold(:,tmp(:)');
            Scur=pinv(KIcur)*X;
            orrc=zeros(3*length(indcur),length(indcur));
            for ic=1:length(indcur),
                [Vy,~,~]=svd(Scur(3*ic-2:3*ic,:));
                orrc(3*ic-2:3*ic,ic)=Vy(:,1);
            end
            Ktildec=KIcur*orrc;
            Stildec=pinv(Ktildec)*X;%Xim;
            Xhatc=Ktildec*Stildec;
            err_resc=norm(X-Xhatc,'fro');
            critsubv(idip)=err_resc/normX+lambda*length(indcur);
        end
        [critsub,ichoix]=min(critsubv);
    end
    
    % choix
    [crit,icrit]=min([crit, critadd, critsub]);
    if icrit==3,
        iic=[1:ichoix-1 ichoix+1:length(Iold)];
        I=Iold(iic);
        tmp=[3*iic-2;3*iic-1;3*iic];
        KI=KIold(:,tmp(:)');
        Jhat=pinv(KI)*X;
        Jblk=zeros(3*length(I),length(I));
        for ic=1:length(I),
            [Vy,~,~]=svd(Jhat(3*ic-2:3*ic,:));
            Jblk(3*ic-2:3*ic,ic)=Vy(:,1);
        end
        ShatI=pinv(KI*Jblk)*X;
        XhatI=KI*Jblk*ShatI; % these are the estimates after eliminating (SBR)
        err_rec=norm(X-XhatI,'fro');
        GOF=1-(err_rec/normX)^2;
    elseif icrit==1,
        I=Iold;
        Jblk=Jblkold;
        ShatI=ShatIold;
        GOF=GOFold;
%        KI=KIold;
%         XhatI=KI*Jblk*ShatI; % these are the estimates at the previous iteration
%         err_reco=norm(X-XhatI,'fro');
%         GOFout=1-(err_reco/normX)^2;
    end
    %end
    
%    err_rec=norm(X-XhatI,'fro');
%    GOF=1-(err_rec/normX)^2;
%     stop1=(GOF>expGOF)&&(icrit~=2);%(GOFo>expGOF)&&(GOFs>expGOF); % all GOF atteint mais pas juste après un ajout ou rejet
%     stop2=(icrit==1);
%     if stop1
%         Iout=I;
%         ShatIout=ShatI;
%         GOFout=GOF;
%     end
%     if stop2 && ~stop1,
%         Jblk=Jblkold;
%     end
%     stop=stop1||stop2;
    stop=(icrit==1);
    %crit=min(critadd,critsub);
    if length(I)>1 && icrit==2
        newlambda=max(0,(err_rec_old-err_rec)/normX);
        lambda=min(lambda,newlambda); % CSBR step (recompute lambda)
    end
    lambda
    %stop=(GOF>expGOF)||(ndipcur>=nbsour);%|(GOF<GOFs);
end
orr=reshape(sum(Jblk,2),3,size(Jblk,1)/3);
% if GOF<GOFs,
% Lf=Lf(:,1:end-3);
% Y1=pinv(Lf)*Xre10;
% Ysol=zeros(3*(ndipcur-1),ndipcur-1);
% for ic=1:ndipcur-1
%     Ysol(3*ic-2:3*ic,ic)=Y1(3*ic-2:3*ic);
% end
% coefs=pinv(Lf*Ysol)*Xreim0;
% %coefs=abs(diag(1./coefs(:,1)))*coefs;
% Xreimh=Lf*Ysol*coefs;
% Xreim=Xreim0-Xreimh;
% GOF=1-(norm(Xreim,'fro')/norm(Xreim0,'fro'))^2
% inddip=inddip(1:end-1);
% end

% indini=1:ndipcur;
% Xreim0r=Xreim0;
% Lft=zeros(M,3,ndipcur);
% Xhc=zeros(M,nhar,ndipcur);
% for idip=indini
%     Lft(:,:,idip)=Lf(:,3*idip-2:3*idip);
%     Xhc(:,:,idip)=Lft(:,:,idip)*Ysol(3*idip-2:3*idip,:);
% end
% Xexp=0;
% inddipper=[];
% while ~isempty(indini)
%     errfin=Inf(1,ndipcur);
%     for idip=indini
%         errfin(idip)=norm(Xreim0r-Xhc(:,:,idip),'fro');
%     end
%     [~,idip]=min(errfin);
%     %idip=indini(1);
%     %errfin
%     %idip=input('idip')
%     inddipper=[inddipper idip];
%     indini=setdiff(indini,idip,'stable');
%     %Xreim0r=Xreim0r-Lft(:,:,idip)*Ysol(3*idip-2:3*idip,idip)*coefs(idip,:);
%     Xexp=Xexp+Xhc(:,:,idip);
%     Xreim0r=Xreim0-Xexp;
%     GOFt=1-(norm(Xreim0-Xexp,'fro')/norm(Xreim0,'fro'))^2;
% end
% inddip=inddip(inddipper);
% Ysols=Ysol;
% tmp=[inddipper*3-2;inddipper*3-1;inddipper*3];
% Ysol=Ysol(tmp(:),:);

