function [idip,Osol,Ssol,GOF]=mySBRR1_AltOpt(V,K,lambda)

lambda=lambda*norm(V,'fro');
ndip=size(K,2)/3;
idip=[];
iset=(1:ndip);
Osol=[];
errk=norm(V,'fro');
errkold=2*errk;
errt=Inf*ones(ndip,1);


while errkold-errk>0 
    k=length(idip);
    id=[3*(idip-1)+1;3*(idip-1)+2;3*idip];
    % ajout d'un bloc regresseur
    for i=iset        
        Kt=[K(:,id(:))*Osol K(:,3*(i-1)+1:3*i)];
        D=pinv(Kt)*V;
        % sol rank 1 with i
        [U,~,~]=svd(D(k+1:k+3,:));
        Osolt=[Osol zeros(3*k,1); zeros(3,k) U(:,1)];
        idt=[id (3*(i-1)+1:3*i)'];
        [Osolt,errt(i)]=OptAlt(V,K(:,idt(:)),Osolt); % une passe d'optimisation alternée
        if (errk-errt(i)-(k+1)*lambda>0) 
            errdiffi=1;
            errti=errt(i);
            while(errdiffi>lambda*1e-1) % itération jusqu'à convergence
                errt(i)=errti;
                [Osolt,errti]=OptAlt(V,K(:,idt(:)),Osolt); % une passe d'optimisation alternée
                errdiffi=errt(i)-errti;
            end
        end
        % ajout pénalisation
        errt(i)=errt(i)+(k+1)*lambda;
        if errt(i)==min(errt)
            imin=i;
            Osolimin=Osolt;
        end
    end
    
     % retrait d'un bloc regresseur
    for j=1:k 
        i=idip(j);
        idipt=setdiff(idip,i,'stable');
        idt=[3*idipt-2;3*idipt-1;3*idipt];
        Osolt=Osol([1:3*(j-1) 3*j+1:3*k],[1:j-1 j+1:k]);        
        [Osolt,errt(i)]=OptAlt(V,K(:,idt(:)),Osolt); % une passe d'optimisation alternée
        if (errt(i)+k*lambda-errk<lambda)
            errdiffi=1+lambda;
            errti=errt(i);
            while(errdiffi>lambda*1e-1) % itération jusqu'à convergence
                errt(i)=errti;
                [Osolt,errti]=OptAlt(V,K(:,idt(:)),Osolt); % une passe d'optimisation alternée
                errdiffi=errt(i)-errti;
            end
        end
        % ajout pénalisation
        errt(i)=errt(i)+(k-1)*lambda;
        if errt(i)<=min(errt)
            imin=i;
            Osolimin=Osolt;
        end
    end
    % choix du meilleur regresseur
    errkold=errk;
    errk=errt(imin);
    if errkold-errk>0
        Osol=Osolimin;
        if ismember(imin,iset)
            idip=[idip imin];
            iset=setdiff(iset,imin);
            %disp(['fin itération - position ajoutée: ', num2str(imin)])
        else
            idip=setdiff(idip,imin,'stable');
            iset=[iset imin];
            %disp(['fin itération - position enlevée: ', num2str(imin)])
        end      
        
    end
    
end

% Un dernier passage d'optim alternée pour la solution finale
errdiffi=lambda;
errti=errkold;
while(errdiffi>lambda*1e-3) % itération jusqu'à convergence
    errt=errti;
    [Osol,errti]=OptAlt(V,K(:,id(:)),Osol); % une passe d'optimisation alternée
    errdiffi=errt-errti;
end


Ksol=K(:,id(:))*Osol;
Ssol=pinv(Ksol)*V;
GOF=1-(norm(V-Ksol*Ssol,'fro')/norm(V,'fro'))^2;

        