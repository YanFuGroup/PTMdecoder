function [FDR,Iid,threshold,finalFDR,p,pp] = ComputeFDR(DecoyType,GroupType,scores,numrst,I,fdrthres)
% compute FDR and return indexes passed the filtering
% input:
%           DecoyType
%               one for each PSM, indicating whether it is a decoy match
%           GroupType
%               one for each PSM, indicating whether it is out of the
%               target group, 0 is wanted
%           scores
%               the score of search engine, sorted
%           numrst
%               the number of results
%           I
%               the index of the PSMs in descend sort
%           fdrthres
%               the FDR threshold set by user
% output:
%           FDR
%               GF (global), SF (separate), TF (transfer).
%               each q-value for each sorted score (the input score)
%           Iid
%               the index of PSM pass the filtering
%           threshold
%               the true score threshold for filtering the PSMs
%           finalFDR
%               the true p-value threshold for filtering the PSMs
%           p
%               the coefficient of fitting in transfer FDR
%           pp
%               the coefficient of lambda fitting

draw_figures = false;   % whether to draw figures of lambda fitting and transfer FDR

II = find(GroupType==0); % The index of the wanted group
numtt = cumsum(DecoyType(II));
numdd = cumsum(~DecoyType(II));
numd = cumsum(~DecoyType);
numt = cumsum(DecoyType);
scores_group = scores(II);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Global FDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDR.GF = numd./numt;
J = (FDR.GF<=fdrthres);
if sum(J)==0
    Iid.GF = [];
    threshold.GF = Inf;
else
    threshold.GF = min(scores(J));
    Iid.GF = (DecoyType  & ~GroupType & scores>=threshold.GF);
end
Iid.GF = I(Iid.GF);
if ~isempty(find(scores==threshold.GF, 1))
    finalFDR.GF = min(FDR.GF(find(scores==threshold.GF)));
else
    finalFDR.GF = 1;
end
% 

[FDR.GF] = Sortfdr(FDR.GF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SeparateFDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(numtt)
    keyboard
end
if numtt(1) == 0
    numtt(1) = 0.000001;
end
Iid.SF = zeros(1,numrst);
FDR.SF = numdd./numtt;

J = (FDR.SF<=fdrthres);
if sum(J)==0
    threshold.SF = Inf;
else
    threshold.SF = min(scores_group(J));
    Iid.SF = ( Iid.SF | ( scores>=threshold.SF) );
end
Iid.SF = Iid.SF & DecoyType & ~GroupType;

Iid.SF = I(Iid.SF);

if ~isempty(find(scores_group==threshold.SF, 1))
    finalFDR.SF = min(FDR.SF(find(scores_group==threshold.SF)));
else
    finalFDR.SF = 1;
end
[FDR.SF] = Sortfdr(FDR.SF);
FDR.SF(FDR.SF>1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TransferredFDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio = numd(II)./numtt;
decoyScores = scores(DecoyType==0);
decoyGroupType = GroupType(DecoyType==0);
III = find(decoyGroupType==0);

if ~isempty(III) && length(III) > 2
    groupDecoyScores = decoyScores(III);
    numGroup = cumsum(~decoyGroupType);
    decoyRatio = numGroup./(1:length(decoyGroupType));
    p = robustfit(groupDecoyScores,decoyRatio(III)); m=p(1); p(1)=p(2); p(2)=m;
    y11 = decoyRatio(III);
    a=decoyRatio(III)>=0;
    groupDecoyScores1=groupDecoyScores(a);
    decoyRatio1=y11(a);
else
     p = [0,0];
end


FDR.TF = ratio.*(scores(II)*p(1)+p(2));


Iid.TF = zeros(1,numrst);
if draw_figures
    if ~isempty(III)
        figure
        plot(groupDecoyScores1,decoyRatio1,'.');
        hold on
        plot(groupDecoyScores1,groupDecoyScores1*p(1)+p(2),'g-','LineWidth',2);
    %     xlim([0,50]);
    else
       figure
       plot(scores,0,'g-','LineWidth',2);
    end
    xlabel('mascot search score','fontsize',12,'fontweight','b');
    ylabel('proportion','fontsize',12,'fontweight','b');
    title('TransferredFDR','fontsize',12,'fontweight','b');
end




for j=1:length(FDR.TF)
    if FDR.TF(j)<0
        FDR.TF(j)=0;
    end
end
J = FDR.TF<=fdrthres;
if sum(J)==0
    threshold.TF = Inf;
else
    threshold.TF = min(scores_group(J));
    Iid.TF = ( Iid.TF | (scores>=threshold.TF) );  
end
Iid.TF = Iid.TF & DecoyType & ~GroupType;
Iid.TF = I(Iid.TF);
if ~isempty(find(scores_group==threshold.TF, 1))
    finalFDR.TF = min(FDR.TF(find(scores_group==threshold.TF)));
else
    finalFDR.TF = 1;
end
[FDR.TF] = Sortfdr(FDR.TF);
FDR.TF(FDR.TF>1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lammada%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lambda fitting plot
ratioo = (numtt-numdd)./(numt(GroupType==0)-numd(GroupType==0));
pp = robustfit(scores_group,ratioo); m=pp(1); pp(1)=pp(2); pp(2)=m;
a = ratioo>=0;
b = scores_group(a);
c = ratioo(a);

if draw_figures
    figure
    plot(b,c,'.');
    hold on
    plot(scores_group,scores_group*pp(1)+pp(2),'g-','LineWidth',2);
    
    xlabel('mascot search score','fontsize',12,'fontweight','b');
    ylabel('proportion','fontsize',12,'fontweight','b');
    title('TransferredFDR','fontsize',12,'fontweight','b');
    % 
    % axes('position',[0.59,0.59,0.3,0.3]);
    % plot(b,c,'.');
    % hold on
    % plot(scores_group,scores_group*pp(1)+pp(2),'g-','LineWidth',2);
    ylim([-0.01,0.01]);
    % 
    % hold off
end

end

%%%% sort fdr %%%%
function [FDR] = Sortfdr(FDR)
for i=length(FDR):-1:2
    if FDR(i-1)>FDR(i)
        FDR(i-1)=FDR(i);
    end
end
end
%%%% sort fdr End %%%%