% Generate non-redundant ions for all IMPs
% Input: 
%   modSites - possible modification positions (numbers) on the peptide sequence
%   massArrangement - arrangement of possible modifications at each position above, each row is a case
%   fixedPosMod - fixed modification
% Output: 
%   vNonRedunTheoryIonMz - the matrix of ion properties determined by all modification combinations, 
%       each row is [m/z, bya type, position index, charge, number of modifications, [by position charge] ion index, ...
%                   0-1 variable indicating whether each IMP can generate this ion]


% Each ion can generate b and y ions
function [vNonRedunTheoryIonMz] = getNonRedunIons(obj,modSites,massArrangement,fixedPosMod)
% obj.m_ionTypes is the array of ion types specified by the user
% 1 is b ion, 11 is b-NH3 ion, 12 is b-H2O ion
% 2 is y ion, 21 is y-NH3 ion, 22 is y-H2O ion
% 3 is a ion, 31 is a-NH3 ion, 32 is a-H2O ion

theoryMz=calculateIonMz(obj,fixedPosMod);% m/z of various by ions without variable modifications
[fragPosNum,maxCharge]=size(theoryMz);
maxCharge=maxCharge/2;

%% Calculate various theoretical ions and their properties
numSite=length(modSites);% number of candidate modification sites
modSites(modSites==0)=1; % N/C terminal modifications are counted on the sequence for convenience
modSites(modSites==length(obj.m_pepSeq)+1)=length(obj.m_pepSeq);
vNonRedunTheoryIonMz=ones(size(massArrangement,1)*(length(obj.m_pepSeq)-1)* ...
    maxCharge*9,6); % information of non-redundant ions
iNRStart=1;% start index for new records in vNonRedunTheoryIonMz
% Localize the modification at the numSite-th site, recalculate the theoretical m/z of the corresponding by ions after adding the modification
for eachArra=1:size(massArrangement,1)% r traverses each modification arrangement
    for i=1:numSite-1% number of modifications on the left of the cleavage site, at least 1, at most number of modifications-1
        b_modMass=sum(massArrangement(eachArra,1:i));
        y_modMass=sum(massArrangement(eachArra,i+1:end));
        nonredundant_b=modSites(i):(modSites(i+1)-1);% position index
        nonredundant_y=(fragPosNum+2-modSites(i+1)):(fragPosNum+1-modSites(i));
        num_nr=length(nonredundant_b);
        if isempty(num_nr)
            % When both N/C terminal and N/C terminal amino acids are modified, cannot distinguish IMPs using ions
            continue;
        end

        for icharge=1:maxCharge% traverse ion charge
            if any(1==obj.m_ionTypes)
                % b_ion=[m/z of site-determining b ions at p charge, same number of 1 (b ion),...
                %       ion index, charge, number of modifications on b side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge+b_modMass)/icharge,8),ones(num_nr,1),...
                    nonredundant_b',icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(11==obj.m_ionTypes)
                % b-NH3 ion, type 11
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),11*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(12==obj.m_ionTypes)
                % b-H2O ion, type 12
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -2*CConstant.hmass-CConstant.omass+b_modMass)/icharge,8),12*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(2==obj.m_ionTypes)
                % y_ion=[m/z of site-determining y ions at p charge, same number of 2 (y ion),...
                %       ion index, charge, number of modifications on y side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge+y_modMass)/icharge,8),2*ones(num_nr,1),...
                    nonredundant_y',icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(21==obj.m_ionTypes)
                % y-NH3 ion, type 21
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge ...
                    -CConstant.nmass-3*CConstant.hmass+y_modMass)/icharge,8),21*ones(num_nr,1),nonredundant_y',...
                    icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(22==obj.m_ionTypes)
                % y-H2O ion, type 22
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge ...
                    -2*CConstant.hmass-CConstant.omass+y_modMass)/icharge,8),22*ones(num_nr,1),nonredundant_y',...
                    icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(3==obj.m_ionTypes)
                % a_ion=[m/z of site-determining a ions at p charge, same number of 3 (a ion),...
                %       ion index, charge, number of modifications on b side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.cmass-CConstant.omass+b_modMass)/icharge,8),3*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(31==obj.m_ionTypes)
                % a-NH3 ion, type 31
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.cmass-CConstant.omass-CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),31*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(32==obj.m_ionTypes)
                % a-H2O ion, type 32
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.cmass-CConstant.omass-CConstant.omass-2*CConstant.hmass+b_modMass)/icharge,8),32*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            %[b+;[b-NH3]+;[b-H2O]+;y+;[y-NH3]+;[y-H2O]+;a+;[a-NH3]+;[a-H2O]+;b++;......]
        end
    end
end
vNonRedunTheoryIonMz(iNRStart:end,:)=[];
vNonRedunTheoryIonMzTemp=vNonRedunTheoryIonMz;% keep the last column for the right half of the matrix
vNonRedunTheoryIonMz=unique(vNonRedunTheoryIonMz(:,1:5),'rows');% remove redundancy

% Add grouping marker
tmpGroupFactor=vNonRedunTheoryIonMz(:,2:4);% only record grouping factors: b,y ion type / fragment ion index / charge
tmpGroupFactor=unique(tmpGroupFactor,'rows','stable');% remove redundancy in original row order to get each feature
[~,IndX]=ismember(vNonRedunTheoryIonMz(:,2:4),tmpGroupFactor,'rows');% IndX is the class number of each row [by, position, charge]
[~,indMassArra]=ismember(vNonRedunTheoryIonMzTemp(:,1:5),vNonRedunTheoryIonMz,'rows');
vNonRedunTheoryIonMz=[vNonRedunTheoryIonMz,zeros(size(vNonRedunTheoryIonMz,1),1+size(massArrangement,1))];
% After column 6, each column represents a IMP, 0-1 variable indicates whether the IMP can generate this ion, 1 means yes
vNonRedunTheoryIonMz((5+vNonRedunTheoryIonMzTemp(:,6))*size(vNonRedunTheoryIonMz,1)+indMassArra)=1;
vNonRedunTheoryIonMz(:,6)=IndX;
end


%% a, b, y ions can all be generated, but only ions with D,E,S,T amino acids can lose water, and only ions with KRNQ amino acids can lose ammonia
% function [vNonRedunTheoryIonMz] = getNonRedunIons(obj,modSites,massArrangement,fixedPosMod)
% % 1 is b ion, 11 is b-NH3 ion, 12 is b-H2O ion
% % 2 is y ion, 21 is y-NH3 ion, 22 is y-H2O ion
% % 3 is a ion, 31 is a-NH3 ion, 32 is a-H2O ion
% 
% theoryMz=calculateIonMz(obj,fixedPosMod);% m/z of b,y ions without variable modifications
% [fragPosNum,maxCharge]=size(theoryMz);
% maxCharge=maxCharge/2;
% [posFirstLossNH3AA,posLastLossNH3AA,posFirstLossH2OAA,posLastLossH2OAA] = ...
%     findNeutralLossStartPosition(obj.m_pepSeq);
% % Theoretical ions and their properties
% numSite=length(modSites);% number of candidate modification sites
% modSites(modSites==0)=1; % N/C terminal modifications are counted on the sequence for convenience
% modSites(modSites==length(obj.m_pepSeq)+1)=length(obj.m_pepSeq);
% vNonRedunTheoryIonMz=ones(size(massArrangement,1)*(length(obj.m_pepSeq)-1)* ...
%     maxCharge*9,6); % information of non-redundant ions
% iNRStart=1;% start index for new records in vNonRedunTheoryIonMz
% % Localize the modification at the numSite-th site, recalculate the theoretical m/z of the corresponding by ions after adding the modification
% for eachArra=1:size(massArrangement,1)% traverse each modification arrangement
%     for i=1:numSite-1% number of modifications on the left of the cleavage site, at least 1, at most number of modifications-1
%         b_modMass=sum(massArrangement(eachArra,1:i));
%         y_modMass=sum(massArrangement(eachArra,i+1:end));
%         nonredundant_b=modSites(i):(modSites(i+1)-1); % position index
%         nonredundant_y=(fragPosNum+2-modSites(i+1)):(fragPosNum+1-modSites(i));
%         num_nr=length(nonredundant_b);
%         if isempty(num_nr)
%             % When both N/C terminal and N/C terminal amino acids are modified, cannot distinguish IMPs using ions
%             continue;
%         end
% 
%         for icharge=1:maxCharge
%             if any(1==obj.m_ionTypes)
%                 % b_ion=[m/z of site-determining b ions at p charge, type indicator 1 (b ion),...
%                 %       ion index, charge, number of modifications on b side, IMP index]
%                 iNREnd=iNRStart+num_nr-1;
%                 vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
%                     [round((theoryMz(nonredundant_b,icharge)*icharge+b_modMass)/icharge,8),ones(num_nr,1),...
%                     nonredundant_b',icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
%                 iNRStart=iNREnd+1;
%             end
%             if any(11==obj.m_ionTypes)
%                 % b-NH3 ion, type 11
%                 for idx=1:num_nr
%                     if nonredundant_b(idx)>=posFirstLossNH3AA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_b(idx),icharge)*icharge ...
%                             -CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),...
%                             11,nonredundant_b(idx),icharge,i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             if any(12==obj.m_ionTypes)
%                 % b-H2O ion, type 12
%                 for idx=1:num_nr
%                     if nonredundant_b(idx)>=posFirstLossH2OAA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_b(idx),icharge)*icharge ...
%                             -2*CConstant.hmass-CConstant.omass+b_modMass)/icharge,8),...
%                             12,nonredundant_b(idx),icharge,i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             if any(2==obj.m_ionTypes)
%                 % y_ion=[m/z of site-determining y ions at p charge, same number of 2 (y ion),...
%                 %       ion index, charge, number of modifications on y side, IMP index]
%                 iNREnd=iNRStart+num_nr-1;
%                 vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
%                     [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge+y_modMass)/icharge,8),2*ones(num_nr,1),...
%                     nonredundant_y',icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
%                 iNRStart=iNREnd+1;
%             end
%             if any(21==obj.m_ionTypes)
%                 % y-NH3 ion, type 21
%                 for idx=1:num_nr
%                     if nonredundant_y(idx)>=length(obj.m_pepSeq)+1-posLastLossNH3AA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_y(idx),icharge+maxCharge)*icharge ...
%                             -CConstant.nmass-3*CConstant.hmass+y_modMass)/icharge,8),...
%                             21,nonredundant_y(idx),icharge,numSite-i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             if any(22==obj.m_ionTypes)
%                 % y-H2O ion, type 22
%                 for idx=1:num_nr
%                     if nonredundant_y(idx)>=length(obj.m_pepSeq)+1-posLastLossH2OAA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_y(idx),icharge+maxCharge)*icharge ...
%                             -2*CConstant.hmass-CConstant.omass+y_modMass)/icharge,8),...
%                             22,nonredundant_y(idx),icharge,numSite-i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             if any(3==obj.m_ionTypes)
%                 % a_ion=[m/z of site-determining a ions at p charge, same number of 3 (a ion),...
%                 %       ion index, charge, number of modifications on b side, IMP index]
%                 iNREnd=iNRStart+num_nr-1;
%                 vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
%                     [round((theoryMz(nonredundant_b,icharge)*icharge ...
%                     -CConstant.cmass-CConstant.omass+b_modMass)/icharge,8),3*ones(num_nr,1),nonredundant_b',...
%                     icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
%                 iNRStart=iNREnd+1;
%             end
%             if any(31==obj.m_ionTypes)
%                 % a-NH3 ion, type 31
%                 for idx=1:num_nr
%                     if nonredundant_b(idx)>=posFirstLossNH3AA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_b(idx),icharge)*icharge ...
%                             -CConstant.cmass-CConstant.omass-CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),...
%                             31,nonredundant_b(idx),icharge,i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             if any(32==obj.m_ionTypes)
%                 % a-H2O ion, type 32
%                 for idx=1:num_nr
%                     if nonredundant_b(idx)>=posFirstLossH2OAA
%                         vNonRedunTheoryIonMz(iNRStart,:)= ...
%                             [round((theoryMz(nonredundant_b(idx),icharge)*icharge ...
%                             -CConstant.cmass-2*CConstant.omass-2*CConstant.hmass+b_modMass)/icharge,8),...
%                             32,nonredundant_b(idx),icharge,i,eachArra];
%                         iNRStart=iNRStart+1;
%                     end
%                 end
%             end
%             % [b+;[b-NH3]+;[b-H2O]+;y+;[y-NH3]+;[y-H2O]+;a+;[a-NH3]+;[a-H2O]+;b++;……],
%             % increase charge in turn, put bya ions and neutral loss ions all in
%         end
%     end
% end
% vNonRedunTheoryIonMz(iNRStart:end,:)=[];
% vNonRedunTheoryIonMzTemp=vNonRedunTheoryIonMz;% keep the last column for the right half of the matrix
% vNonRedunTheoryIonMz=unique(vNonRedunTheoryIonMz(:,1:5),'rows');% remove redundancy
% 
% %% Add grouping marker
% tmpGroupFactor=vNonRedunTheoryIonMz(:,2:4);% only record grouping factors: by ion type / fragment ion index / charge
% tmpGroupFactor=unique(tmpGroupFactor,'rows','stable');% remove redundancy in original row order to get each feature
% [~,IndX]=ismember(vNonRedunTheoryIonMz(:,2:4),tmpGroupFactor,'rows');% IndX is the class number of each row [by, position, charge]
% [~,indMassArra]=ismember(vNonRedunTheoryIonMzTemp(:,1:5),vNonRedunTheoryIonMz,'rows');
% vNonRedunTheoryIonMz=[vNonRedunTheoryIonMz,zeros(size(vNonRedunTheoryIonMz,1),1+size(massArrangement,1))];
% % After column 6, each column represents a IMP, 0-1 variable indicates whether this IMP can generate this ion, 1 means yes
% vNonRedunTheoryIonMz((5+vNonRedunTheoryIonMzTemp(:,6))*size(vNonRedunTheoryIonMz,1)+indMassArra)=1;% directly find the element position sorted by column
% vNonRedunTheoryIonMz(:,6)=IndX;
% end
% 
% function [posFirstLossNH3AA,posLastLossNH3AA,posFirstLossH2OAA,posLastLossH2OAA] = ...
%     findNeutralLossStartPosition(pepSeq)
% % Find the first position on N, C terminal that can have neutral loss
% % Input: pepSeq is the peptide sequence
% % Output: posFirstLossNH3AA - the first position that can lose ammonia
% %         posLastLossNH3AA - the last position that can lose ammonia
% %         posFirstLossH2OAA - the first position that can lose water
% %         posLastLossH2OAA - the last position that can lose water
% %         If the corresponding amino acid is not found, first is length(pepSeq)+1, last is 0
% posFirstLossNH3AA=length(pepSeq)+1;
% posLastLossNH3AA=0;
% posFirstLossH2OAA=length(pepSeq)+1;
% posLastLossH2OAA=0;
% 
% lossNH3AA={'K','R','N','Q'};
% lossH2OAA={'D','E','S','T'};
% 
% % Update the starting position of ammonia-losing amino acids on both sides
% for idx=1:length(lossNH3AA)
%     findRes=strfind(pepSeq,lossNH3AA{idx});
%     if ~isempty(findRes)
%         if findRes(1)<posFirstLossNH3AA
%             posFirstLossNH3AA=findRes(1);
%         end
%         if findRes(end)>posLastLossNH3AA
%             posLastLossNH3AA=findRes(end);
%         end
%     end
% end
% % Update the starting position of water-losing amino acids on both sides
% for idx=1:length(lossH2OAA)
%     findRes=strfind(pepSeq,lossH2OAA{idx});
%     if ~isempty(findRes)
%         if findRes(1)<posFirstLossH2OAA
%             posFirstLossH2OAA=findRes(1);
%         end
%         if findRes(end)>posLastLossH2OAA
%             posLastLossH2OAA=findRes(end);
%         end
%     end
% end
% end