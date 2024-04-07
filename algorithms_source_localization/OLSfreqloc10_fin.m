function [I,orr,ShatI,GOF] = OLSfreqloc10_fin(X,K,expGOF)
% SBRFREQLOC8 implements a block-SBR regression under rank 1 constranints
% V=data matrix, one signal per row
% L=lead-field
% expGOF=expected GOF
% indloc=source indices (positions)
% orr=source orientations
% Stilde=source amplitudes
% GOF=GOF

if expGOF>1
    nbsour=expGOF;
else
    nbsour=Inf;
end

normX=norm(X,'fro');

crit=Inf;
GOF=0;
ndip=size(K,2)/3;

Xhat=[];
I=[];

stop=0;

KI=[];
Jblk=[];

while ~stop
    ndipcur=length(I)+1;
    err=Inf(1,ndip);
    err2=Inf(1,ndip);
    I2=I;
    Jblk2=Jblk;
    KI2=KI;
    for idip=setdiff([1:ndip],I)%
        ki=K(:,3*idip-2:3*idip);
        KIt=[KI ki];
        JIt=pinv(KIt)*X;
        Jblkt=zeros(3*ndipcur,ndipcur);
        for ic=1:ndipcur,
            [U,~,~]=svd(JIt(3*ic-2:3*ic,:));
            Jblkt(3*ic-2:3*ic,ic)=U(:,1);
        end
        ShatIt=pinv(KIt*Jblkt)*X;
        Xhatit=KIt*Jblkt*ShatIt;
        err(idip)=norm(X-Xhatit,'fro');
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
    XhatI=KI*Jblk*ShatI;
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
    if errm2<errm,
        I=I2;
        Jblk=Jblk2;
        ShatI=ShatI2;
        KI=KI2;
        XhatI=XhatI2;
    end
    err_rec=norm(X-XhatI,'fro');
    GOF=1-(err_rec/normX)^2
  
    stop=(GOF>expGOF);
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

