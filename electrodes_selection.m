function [velecf1,velecf0,ixf1, ixf0, LFself1,LFself0] = electrodes_selection(LF,velecf,indf1,indf0)
%ELECTRODES_SELECTION Select only electrodes with resolute peaks 
% (considering neighboring information) at the frequencies of interest 
% f1 = 1.2,2.4,3.6,4.8 Hz and f0 = 6,12,18,24 Hz

%% Electrodes with peaks at frequencies of face individuation response 1.2,2.4,3.6,4.8 Hz
q=1;
velecf1=[];
ixf1=[];
LFself1=[];
for j=1:size(velecf,1)
    if  (abs(velecf(j,indf1(1)-1)) < abs(velecf(j,indf1(1)))) && (abs(velecf(j,indf1(1))) > abs(velecf(j,indf1(1)+1))) ...
        || (abs(velecf(j,indf1(2)-1)) < abs(velecf(j,indf1(2)))) && (abs(velecf(j,indf1(2))) > abs(velecf(j,indf1(2)+1))) ...
        || (abs(velecf(j,indf1(3)-1)) < abs(velecf(j,indf1(3)))) && (abs(velecf(j,indf1(3))) > abs(velecf(j,indf1(3)+1)))
        velecf1(q,:)=velecf(j,:);
        ixf1(q,1)=j;
        q=q+1;
    end
end
LFself1=LF(ixf1,:);

%% Electrodes with peaks at frequencies of general visual response 6,12,18,24 Hz
q=1;
velecf0=[];
ixf0=[];
LFself0=[];
for j=1:size(velecf,1)
    if  (abs(velecf(j,indf0(1)-1)) < abs(velecf(j,indf0(1)))) && (abs(velecf(j,indf0(1))) > abs(velecf(j,indf0(1)+1))) ...
        || (abs(velecf(j,indf0(2)-1)) < abs(velecf(j,indf0(2)))) && (abs(velecf(j,indf0(2))) > abs(velecf(j,indf0(2)+1))) ...
        || (abs(velecf(j,indf0(3)-1)) < abs(velecf(j,indf0(3)))) && (abs(velecf(j,indf0(3))) > abs(velecf(j,indf0(3)+1)))
        velecf0(q,:)=velecf(j,:);
        ixf0(q,1)=j;
        q=q+1;
    end
end
LFself0=LF(ixf0,:);
%LFself0=LF_Norm(ixf0,:);
end

