function [indloc,orr,Stilde,GOF]=allmusic(V,L,expGOF,type,datared,critere)
% ALLMUSIC implements different flavours of MUSIC, given by type
% V=data matrix, one signal per row
% L=lead-field
% expGOF=expected GOF
% type='plain','rec','rap','trap', default 'plain'
% datared= 1 or 0 datareduction (source/noise subspaces) default 0
% critere= data reduction method: 'aic': Akaike, 'mdl':min desc length,
%          numeric<1: expGOF, numeric>1: given signal space dimension; if
%          not given (default): keep all eigenvalues>0
% indloc=source indices (positions)
% orr=source orientations
% Stilde=source amplitudes
% GOF=GOF
% inspired from Moiseev 2011, 2021

if nargin<4,
    type='plain';datared=0;
elseif nargin<5
    datared=0;
end

[U,S,~]=svd(V); % signal space

if datared==1,
    if nargin<6,
        [r,~] = aicmdl(V);
    else
        if isnumeric(critere)
            if critere>=1,
                r=critere;
                %expGOF=sum(diag(S(1:r,1:r)).^2)/sum(diag(S).^2);
            else
                expGOF=critere;
                r=find(cumsum(diag(S).^2)/sum(diag(S).^2)>expGOF,1,'first');
                %expGOF=1;
            end
        else
            [r,~] = aicmdl(V,critere);
        end
    end
    %Usave=U;
    %overestimate r
    %r=ceil(1.2*r);
    r=min(r,min(size(V)));
    U=U(:,1:r);
    Usave=U;
end
[m,p]=size(L);
p=p/3;
if p~=round(p)
    error('incorrect leadfield - groups of 3 columns')
end
GOF=0;
if strcmp(type,'plain'),
    % plain music
    P=U*U';
    for iloc=1:p,
        Lcur=L(:,iloc*3-2:iloc*3);
        [o,mut]=eig(Lcur'*P*Lcur,Lcur'*Lcur);
        o=real(o);
        mud=abs(diag(mut));
        [mu(iloc),imud]=max(mud);
        orr(:,iloc)=o(:,imud)/norm(o(:,imud));
    end
    
    [~,inds]=sort(mu,'descend');
    Lfin=[];
    isour=1;
    
    while (GOF<expGOF)%&&(isour<=r)
        Lfin=[Lfin L(:,inds(isour)*3-2:inds(isour)*3)*orr(:,inds(isour))];
        Stilde=pinv(Lfin)*V;
        Vhat=Lfin*Stilde;
        GOF=1-(norm(V-Vhat,'fro')/norm(V,'fro'))^2;
        isour=isour+1;
    end
    indloc=(inds(1:isour-1));
    orr=orr(:,indloc);
elseif strcmp(type,'rap'),
    Ltmp=[];
    isour=1;
    indloc=[];
    orr=[];
    ind1=[];
    while (GOF<expGOF)&&(cond(Ltmp'*Ltmp)<10^20)&&(isour<=r)
        if isour==1
            Q=eye(m);
        else
            Q=eye(m)-Ltmp*(Ltmp'*Ltmp)^(-1)*Ltmp';
        end
        Sk=Q*U;
        P=Sk*(Sk'*Sk)^(-1)*Sk';
        mu=zeros(1,p);
        for iloc=setdiff(1:p,ind1),
            Lcur=L(:,iloc*3-2:iloc*3);
            [o,mut]=eig((Q*Lcur)'*P*(Q*Lcur),(Q*Lcur)'*(Q*Lcur));
            o=real(o);
%             mu(iloc)=abs(mut(1,1));
%             orrt(:,iloc)=o(:,1)/norm(o(:,1));
            mud=abs(diag(mut));
            [mu(iloc),imud]=max(mud);
            orrt(:,iloc)=o(:,imud)/norm(o(:,imud));
        end
        [Mmu,ind1]=max(mu);
        Ltmp=[Ltmp L(:,3*ind1-2:3*ind1)*orrt(:,ind1)];
        if cond(Ltmp'*Ltmp)<10^20
            isour=isour+1;
            Stilde=pinv(Ltmp)*V;
            Vhat=Ltmp*Stilde;
            indloc=[indloc,ind1];
            orr=[orr orrt(:,ind1)];
        end
        GOF=1-(norm(V-Vhat,'fro')/norm(V,'fro'))^2;  
    end
elseif strcmp(type,'rec'),
    Ltmp=[];
    isour=1;
    indloc=[];
    orr=[];
    while (GOF<expGOF)&&(cond(Ltmp'*Ltmp)<10^10)&&(isour<=r)
        if isour==1
            Q=eye(m);
        else
            Q=eye(m)-Ltmp*(Ltmp'*Ltmp)^(-1)*Ltmp';
        end
        P=U*U';
        for iloc=1:p,
            Lcur=L(:,iloc*3-2:iloc*3);
            [o,mut]=eig((Q*Lcur)'*P*(Q*Lcur),(Q*Lcur)'*(Q*Lcur));
            o=real(o);
%             mu(iloc)=abs(mut(1,1));
%             orrt(:,iloc)=o(:,1)/norm(o(:,1));
            mud=abs(diag(mut));
            [mu(iloc),imud]=max(mud);
            orrt(:,iloc)=o(:,imud)/norm(o(:,imud));
        end
        [Mmu,ind1]=max(mu);
        isour=isour+1;
        Ltmp=[Ltmp L(:,3*ind1-2:3*ind1)*orrt(:,ind1)];
        Stilde=pinv(Ltmp)*V;
        Vhat=Ltmp*Stilde;
        GOF=1-(norm(V-Vhat,'fro')/norm(V,'fro'))^2;
        indloc=[indloc,ind1];
        orr=[orr orrt(:,ind1)];
    end
elseif strcmp(type,'trap'),
    extrar=0;
    Ltmp=[];
    isour=1;
    indloc=[];
    orr=[];
    ind1=[];
    while (GOF<expGOF)&&(cond(Ltmp'*Ltmp)<10^20)&&(isour<=r)
        if isour==1
            Q=eye(m);
        else
            Q=eye(m)-Ltmp*(Ltmp'*Ltmp)^(-1)*Ltmp';
        end
        Sk=Q*U(:,1:r+extrar+1-isour);
        P=Sk*(Sk'*Sk)^(-1)*Sk';
        mu=zeros(1,p);
        for iloc=setdiff(1:p,ind1),
            Lcur=L(:,iloc*3-2:iloc*3);
            [o,mut]=eig((Q*Lcur)'*P*(Q*Lcur),(Q*Lcur)'*(Q*Lcur));
            %[o,mut]=eig(Lcur'*(Q'*P*Q)*Lcur,Lcur'*(Q'*Q)*Lcur);
            o=real(o);
%             mu(iloc)=abs(mut(1,1));
%             orrt(:,iloc)=o(:,1)/norm(o(:,1));
            mud=abs(diag(mut));
            [mu(iloc),imud]=max(mud);
            orrt(:,iloc)=o(:,imud)/norm(o(:,imud));
        end
        [Mmu,ind1]=max(mu);
        Ltmp=[Ltmp L(:,3*ind1-2:3*ind1)*orrt(:,ind1)];
        if cond(Ltmp'*Ltmp)<10^20
            isour=isour+1;
            Stilde=pinv(Ltmp)*V;
            Vhat=Ltmp*Stilde;
            indloc=[indloc,ind1];
            orr=[orr orrt(:,ind1)];
        end
        GOF=1-(norm(V-Vhat,'fro')/norm(V,'fro'))^2; 
    end
    
end


