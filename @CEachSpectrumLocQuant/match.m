% Iterate through various modifications to find the union of experimental spectrum peak sets that can match the mass-to-charge ratio
% Input: 
%   expPeaks (N x 2 double)
%       experimental spectrum peaks [m/z, intensity]
%   vNonRedunTheoryIonMz (L x 1 double or L x M double)
%       ion info matrix, first column is m/z
% Output: 
%   vMatchedExpPeaks (K x 2 double)
%       matched [ion_index, intensity]
%       ion_index is the row index in vNonRedunTheoryIonMz

% Direct matching, return [matched vNonRedunTheoryIonMz number intensity] pair
function [vMatchedExpPeaks]=match(obj,expPeaks,vNonRedunTheoryIonMz)
vMatchedExpPeaks=[(1:size(vNonRedunTheoryIonMz,1))',zeros(size(vNonRedunTheoryIonMz,1),1)];
for idxExpPeak=1:size(expPeaks,1)
    iMatchedNonRedun=find(abs(vNonRedunTheoryIonMz(:,1)-expPeaks(idxExpPeak,1)) ...
        <=obj.m_ms2_tolerance);%Row number matched in vNonRedunTheoryIonMz
    for idxMatched=1:length(iMatchedNonRedun)
        vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2) = ...
            vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2)+expPeaks(idxExpPeak,2);
    end
end
vMatchedExpPeaks(vMatchedExpPeaks(:,2)==0,:)=[];
end

% Add dehydration and deamination to the intensity of spectrum peaks without neutral loss, return [matched vNonRedunTheoryIonMz number intensity] pair
% function [vMatchedExpPeaks]=match(~,expPeaks,vNonRedunTheoryIonMz)
% vMatchedExpPeaks=[(1:size(vNonRedunTheoryIonMz,1))',zeros(size(vNonRedunTheoryIonMz,1),1)];
% ionsWithoutNeutralLoss=zeros(size(expPeaks,1),2);%Complete the spectrum peaks with neutral loss, can only match bya
% idxWithoutLoss=1;
% for idxExpPeak=1:size(expPeaks,1)
%     iMatchedNonRedun=find(abs(vNonRedunTheoryIonMz(:,1)-expPeaks(idxExpPeak,1)) ...
%         <=obj.m_ms2_tolerance);%Row number matched in vNonRedunTheoryIonMz
%     for idxMatched=1:length(iMatchedNonRedun)
%         if any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[11,21,31])%NH3
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (CConstant.nmass+3*CConstant.hmass)/vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         elseif any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[12,22,32])%H2O
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (2*CConstant.hmass+CConstant.omass)/vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         else
%             vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2) = ...
%                 vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2)+expPeaks(idxExpPeak,2);
%         end        
%     end
% end
% 
% for idxExpPeak=1:idxWithoutLoss-1%Iterate through the experimental spectrum ionsWithoutNeutralLoss with completed neutral loss
%     iMatchedNonRedun=find(abs(vNonRedunTheoryIonMz(:,1)-ionsWithoutNeutralLoss(idxExpPeak,1)) ...
%         <=obj.m_ms2_tolerance);
%     for idxMatched=1:length(iMatchedNonRedun)
%         if any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[1,2,3])
%             vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2) = ...
%                 vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2)+ionsWithoutNeutralLoss(idxExpPeak,2);
%         end
%     end
% end
% vMatchedExpPeaks(vMatchedExpPeaks(:,2)==0,:)=[];
% end

% Add dehydration, deamination and a-ion to the intensity of by without neutral loss, return [matched vNonRedunTheoryIonMz number intensity] pair
% function [vMatchedExpPeaks]=match(~,expPeaks,vNonRedunTheoryIonMz)
% vMatchedExpPeaks=[(1:size(vNonRedunTheoryIonMz,1))',zeros(size(vNonRedunTheoryIonMz,1),1)];
% ionsWithoutNeutralLoss=zeros(size(expPeaks,1),2);%Complete the spectrum peaks with neutral loss, can only match by
% idxWithoutLoss=1;
% for idxExpPeak=1:size(expPeaks,1)
%     iMatchedNonRedun=find(abs(vNonRedunTheoryIonMz(:,1)-expPeaks(idxExpPeak,1)) ...
%         <=obj.m_ms2_tolerance);%Row number matched in vNonRedunTheoryIonMz
%     for idxMatched=1:length(iMatchedNonRedun)
%         if any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[11,21])%NH3
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (CConstant.nmass+3*CConstant.hmass)/vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         elseif any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[12,22])%H2O
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (2*CConstant.hmass+CConstant.omass)/vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         elseif isequal(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2),3)%a-ion
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (CConstant.cmass+CConstant.omass)/vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         elseif isequal(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2),31)%a-ion with NH3 loss
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (CConstant.cmass+CConstant.omass+CConstant.nmass+3*CConstant.hmass)/...
%                 vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         elseif isequal(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2),32)%a-ion with H2O loss
%             ionsWithoutNeutralLoss(idxWithoutLoss,:) = [expPeaks(idxExpPeak,1)+...
%                 (CConstant.cmass+CConstant.omass+2*CConstant.hmass+CConstant.omass)/...
%                 vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),4),...
%                 expPeaks(idxExpPeak,2)];
%             idxWithoutLoss = idxWithoutLoss+1;
%         else
%             vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2) = ...
%                 vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2)+expPeaks(idxExpPeak,2);
%         end        
%     end
% end
% 
% for idxExpPeak=1:idxWithoutLoss-1%Iterate through the experimental spectrum ionsWithoutNeutralLoss with completed neutral loss
%     iMatchedNonRedun=find(abs(vNonRedunTheoryIonMz(:,1)-ionsWithoutNeutralLoss(idxExpPeak,1)) ...
%         <=obj.m_ms2_tolerance);
%     for idxMatched=1:length(iMatchedNonRedun)
%         if any(vNonRedunTheoryIonMz(iMatchedNonRedun(idxMatched),2)==[1,2])
%             vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2) = ...
%                 vMatchedExpPeaks(iMatchedNonRedun(idxMatched),2)+ionsWithoutNeutralLoss(idxExpPeak,2);
%         end
%     end
% end
% 
% vMatchedExpPeaks(vMatchedExpPeaks(:,2)==0,:)=[];
% end