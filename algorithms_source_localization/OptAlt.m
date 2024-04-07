function [O,err]=OptAlt(V,K,O)
% passe d'optimisation alternée

k=size(K,2)/3;
id=(1:k);
id3=(1:3*k);

for i=1:k
    idsi=setdiff(id,i);
    id3i=(3*(i-1)+1:3*i);
    id3si=setdiff(id3,id3i);
    Kt=[K(:,id3si)*O(id3si,idsi) K(:,3*(i-1)+1:3*i)];
    D=pinv(Kt)*V;
    % sol rank 1 with i
    [U,~,~]=svd(D(k:k+2,:));
    O(id3i,i)=U(:,1);
end

err=norm(V-(K*O)*pinv(K*O)*V,'fro');

