function [vNonRedunTheoryIonMz] = getNonRedunIons(ctx,modSites,massArrangement,fixedPosMod)
% Generate non-redundant ions for all IMPs
% Inputs:
%   ctx (struct)
%       Required fields: m_pepSeq, m_ionTypes, m_iCharge.
%       m_ionTypes is a vector of user-configured ion types, such as [1,2,3,11,12,...], where 
%       1 is b ion, 11 is b-NH3 ion, 12 is b-H2O ion
%       2 is y ion, 21 is y-NH3 ion, 22 is y-H2O ion
%       3 is a ion, 31 is a-NH3 ion, 32 is a-H2O ion
%   modSites (1 x S double/int)
%       Candidate modification-site positions.
%   massArrangement (M x S double)
%       Candidate mass arrangements for each peptidoform.
%   fixedPosMod (K x 3 cell)
%       Fixed modification list [position, mod_name, mod_mass].
% Outputs:
%   vNonRedunTheoryIonMz (L x T double)
%       Non-redundant ion table. Columns:
%       [m/z, ion_type, ion_pos, charge, mod_count_side, ion_group_idx, IMP-membership(0-1 variable indicating whether the ion belongs to the IMP)].

theoryMz=CMS2MassCalculator.calculateIonMz(ctx,fixedPosMod);
[fragPosNum,maxCharge]=size(theoryMz);
maxCharge=maxCharge/2;

numSite=length(modSites);
modSites(modSites==0)=1;
modSites(modSites==length(ctx.m_pepSeq)+1)=length(ctx.m_pepSeq);
vNonRedunTheoryIonMz=ones(size(massArrangement,1)*(length(ctx.m_pepSeq)-1)* ...
    maxCharge*9,6);
iNRStart=1;

for eachArra=1:size(massArrangement,1)
    for i=1:numSite-1
        b_modMass=sum(massArrangement(eachArra,1:i));
        y_modMass=sum(massArrangement(eachArra,i+1:end));
        nonredundant_b=modSites(i):(modSites(i+1)-1);
        nonredundant_y=(fragPosNum+2-modSites(i+1)):(fragPosNum+1-modSites(i));
        num_nr=length(nonredundant_b);
        if isempty(num_nr)
            continue;
        end

        for icharge=1:maxCharge
            % Build selected ion channels according to user-configured ion types.
            if any(1==ctx.m_ionTypes)
                % b_ion=[m/z of site-determining b ions at p charge, same number of 1 (b ion),...
                %       ion index, charge, number of modifications on b side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge+b_modMass)/icharge,8),ones(num_nr,1),...
                    nonredundant_b',icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(11==ctx.m_ionTypes)
                % b-NH3 ion, type 11
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),11*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(12==ctx.m_ionTypes)
                % b-H2O ion, type 12
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -2*CConstant.hmass-CConstant.omass+b_modMass)/icharge,8),12*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(2==ctx.m_ionTypes)
                % y_ion=[m/z of site-determining y ions at p charge, same number of 2 (y ion),...
                %       ion index, charge, number of modifications on y side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge+y_modMass)/icharge,8),2*ones(num_nr,1),...
                    nonredundant_y',icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(21==ctx.m_ionTypes)
                % y-NH3 ion, type 21
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge ...
                    -CConstant.nmass-3*CConstant.hmass+y_modMass)/icharge,8),21*ones(num_nr,1),nonredundant_y',...
                    icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(22==ctx.m_ionTypes)
                % y-H2O ion, type 22
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_y,(icharge+maxCharge))*icharge ...
                    -2*CConstant.hmass-CConstant.omass+y_modMass)/icharge,8),22*ones(num_nr,1),nonredundant_y',...
                    icharge*ones(num_nr,1),(numSite-i)*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(3==ctx.m_ionTypes)
                % a_ion=[m/z of site-determining a ions at p charge, same number of 3 (a ion),...
                %       ion index, charge, number of modifications on b side, IMP index]
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.cmass-CConstant.omass+b_modMass)/icharge,8),3*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(31==ctx.m_ionTypes)
                % a-NH3 ion, type 31
                iNREnd=iNRStart+num_nr-1;
                vNonRedunTheoryIonMz(iNRStart:iNREnd,:)= ...
                    [round((theoryMz(nonredundant_b,icharge)*icharge ...
                    -CConstant.cmass-CConstant.omass-CConstant.nmass-3*CConstant.hmass+b_modMass)/icharge,8),31*ones(num_nr,1),nonredundant_b',...
                    icharge*ones(num_nr,1),i*ones(num_nr,1),eachArra*ones(num_nr,1)];
                iNRStart=iNREnd+1;
            end
            if any(32==ctx.m_ionTypes)
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

% Deduplicate ion rows by [m/z, type, pos, charge, mod_count], then append:
% - group index for (type,pos,charge)
% - IMP membership indicator columns
vNonRedunTheoryIonMzTemp=vNonRedunTheoryIonMz;
vNonRedunTheoryIonMz=unique(vNonRedunTheoryIonMz(:,1:5),'rows');

tmpGroupFactor=vNonRedunTheoryIonMz(:,2:4);
tmpGroupFactor=unique(tmpGroupFactor,'rows','stable');
[~,IndX]=ismember(vNonRedunTheoryIonMz(:,2:4),tmpGroupFactor,'rows');% IndX is the class number of each row [by, position, charge]
[~,indMassArra]=ismember(vNonRedunTheoryIonMzTemp(:,1:5),vNonRedunTheoryIonMz,'rows');
vNonRedunTheoryIonMz=[vNonRedunTheoryIonMz,zeros(size(vNonRedunTheoryIonMz,1),1+size(massArrangement,1))];
% After column 6, each column represents a IMP, 0-1 variable indicates whether the IMP can generate this ion, 1 means yes
vNonRedunTheoryIonMz((5+vNonRedunTheoryIonMzTemp(:,6))*size(vNonRedunTheoryIonMz,1)+indMassArra)=1;
vNonRedunTheoryIonMz(:,6)=IndX;
end
