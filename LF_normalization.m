function LF_Norm = LF_normalization(LF)
%LF_NORMALIZATION Normalize the LF data to keep the same scale on them 

[nel,ndip3] = size(LF);
ndip = ndip3/3;
for idip = 1:ndip
    tmp = LF(:,3*idip-2:3*idip);
    tmpCov = tmp'*tmp;
    LF_Norm(:,3*idip-2:3*idip) = tmp*tmpCov^(-1/2);
end

end

