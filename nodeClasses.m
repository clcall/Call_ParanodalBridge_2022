function ALL = nodeClasses(ALL)
    cond = {'ctrl','cupr'};
    ALL.stats.nodeClassLasts = [];
    ALL.stats.nodeClassBslns = [];
    for c = 1:2
        cells = fieldnames(ALL.(cond{c}));
        % INITITALIZE OUTPUTS
        ALL.stats.NodeClassLnthAvg_bsln.(cond{c}) = [];
        ALL.stats.NodeClassLnthAvg_last.(cond{c}) = [];
        ALL.stats.NodeClassTots_bsln.(cond{c}) = [];
        ALL.stats.NodeClassTots_last.(cond{c}) = [];
        ALL.stats.lostIntClass.(cond{c}) = [];
        ALL.stats.lostIntClassProp.(cond{c}) = [];
        ALL.stats.i2tc_LnthChng.(cond{c}) = [];
        ALL.stats.i2t_LnthChngs.(cond{c}) = [];
        ALL.stats.i2c_LnthChngs.(cond{c}) = [];
        ALL.stats.i2i_LnthChngs.(cond{c}) = [];
        ALL.stats.t2c_LnthChngs.(cond{c}) = [];
        ALL.stats.isolDynam.(cond{c}) = [];
        ALL.stats.isolStbl.(cond{c}) = [];
        % Parse data cell by cell
        for k = 1:length(cells)
            nodes = ALL.(cond{c}).(cells{k}).nodeClassPerTP;
            lnths = ALL.(cond{c}).(cells{k}).sheathLengthPerTP;
            temp = lnths==0;
            lostIdx = any(temp,2);
            mostnbrs = max(nodes(lostIdx,:),[],2,'omitnan');
            numIsol = sum(max(nodes,[],2)==3);
            numTerm = sum(max(nodes,[],2)==4);
            numCont = sum(max(nodes,[],2)==5);
            isolClass = sum(any(mostnbrs==3,2));
            termClass = sum(any(mostnbrs==4,2));
            contClass = sum(any(mostnbrs==5,2));
            % numbers of internodes lost per class
            ALL.stats.lostIntClass.(cond{c}) = [ALL.stats.lostIntClass.(cond{c}); isolClass termClass contClass];
            % proportions of internodes per class
            isolProp = isolClass./numIsol;
            termProp = termClass./numTerm;
            contProp = contClass./numCont;
            ALL.stats.lostIntClassProp.(cond{c}) = [ALL.stats.lostIntClassProp.(cond{c}); isolProp termProp contProp];
            % parsing stable and dynamic internodes
            i2i = max(nodes,[],2)==3;
            i2i_lnths = lnths(i2i,:);
            i2i_netchng = i2i_lnths(:,end) - i2i_lnths(:,1);
            
            i2t = nodes(:,1)==3 & nodes(:,end)==4;
            i2t_lnths = lnths(i2t,:);
            i2t_netchng = i2t_lnths(:,end) - i2t_lnths(:,1);
            
            i2c = nodes(:,1)==3 & nodes(:,end)==5;
            i2c_lnths = lnths(i2c,:);
            i2c_netchng = i2c_lnths(:,end) - i2c_lnths(:,1);
            
            t2t = nodes(:,1)==4 & max(nodes,[],2)==4;
            t2t_lnths = lnths(t2t,:);
            t2t_netchng = t2t_lnths(:,end) - t2t_lnths(:,1);
            
            t2c = nodes(:,1)==4 & nodes(:,end)==5;
            t2c_lnths = lnths(t2c,:);
            t2c_netchng = t2c_lnths(:,end) - t2c_lnths(:,1);
            
            c2c = nodes(:,1)==5 & max(nodes,[],2)==5;
            c2c_lnths = lnths(c2c,:);
            c2c_netchng = c2c_lnths(:,end) - c2c_lnths(:,1);
            
            % assigning length changes to ALL struct
            ALL.stats.i2t_LnthChngs.(cond{c}) = [ALL.stats.i2t_LnthChngs.(cond{c}); i2t_netchng];
            ALL.stats.i2c_LnthChngs.(cond{c}) = [ALL.stats.i2c_LnthChngs.(cond{c}); i2c_netchng];
            ALL.stats.t2c_LnthChngs.(cond{c}) = [ALL.stats.t2c_LnthChngs.(cond{c}); t2c_netchng];
            ALL.stats.isolDynam.(cond{c}) = [ALL.stats.isolDynam.(cond{c}); sum(i2t_netchng>0) + sum(i2c_netchng>0),sum(i2t_netchng<0) + sum(i2c_netchng<0)];
            ALL.stats.isolStbl.(cond{c}) = [ALL.stats.isolStbl.(cond{c}); sum(i2i_netchng>0),sum(i2i_netchng<0)];
            ALL.stats.i2tc_LnthChng.(cond{c}) = [ALL.stats.i2tc_LnthChng.(cond{c}); mean(i2t_netchng) mean(i2c_netchng)];
            ALL.stats.i2i_LnthChngs.(cond{c}) = [ALL.stats.i2i_LnthChngs.(cond{c}); i2i_netchng];
            % get lengths and total numbers of internodes for each
            %       designation at d00 and d14
            days = fieldnames(ALL.(cond{c}).(cells{k}).raw);
            types = [];
            lnths = [];
            for d = 1:length(days)
                nodes = ALL.(cond{c}).(cells{k}).raw.(days{d}).Sheath_Info.node_class;
                lnthsTemp = ALL.(cond{c}).(cells{k}).raw.(days{d}).Sheath_Info.intLengths;
                temp = [sum(nodes==1);...
                        sum(nodes==2);...
                        sum(nodes==3);...
                        sum(nodes==4);...
                        sum(nodes==5)];
                types = [types temp];
                lnths = [lnths lnthsTemp];
                ALL.(cond{c}).(cells{k}).NodeClassTots = types;
                % calculate average lengths per class of internode at d00 and d14
                if d==1
                     ALL.stats.NodeClassLnthAvg_bsln.(cond{c}) = [ALL.stats.NodeClassLnthAvg_bsln.(cond{c});...
                                                       mean(i2i_lnths(:,1)) mean(t2t_lnths(:,1)) mean(c2c_lnths(:,1)) mean(i2t_lnths(:,1)) mean(i2c_lnths(:,1)) mean(t2c_lnths(:,1))];
                elseif d==length(days)
                    if contains(days{end},'d14') || contains(days{end},'d15')
                        ALL.stats.NodeClassLnthAvg_last.(cond{c}) = [ALL.stats.NodeClassLnthAvg_last.(cond{c});...
                                                       mean(i2i_lnths(:,end)) mean(t2t_lnths(:,end)) mean(c2c_lnths(:,end)) mean(i2t_lnths(:,end)) mean(i2c_lnths(:,end)) mean(t2c_lnths(:,end))];
                    else
                        ALL.stats.NodeClassLnthAvg_last.(cond{c}) = [ALL.stats.NodeClassLnthAvg_last.(cond{c});...
                                                       NaN(1,6)];
                    end
                end
            end
        end
    end
    %% average lengths at d00 and d14
    Avgs = [mean(ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,1:3),'omitnan'); mean(ALL.stats.NodeClassLnthAvg_last.ctrl(:,1:3),'omitnan')]; 
    Sems = [calcSEM(ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,1:3)); calcSEM(ALL.stats.NodeClassLnthAvg_last.ctrl(:,1:3))]; 

    colors = {[255 247 147]...
              [215 71 90]...
              [74 254 255]};
    for i=1:3
        colors{i} = colors{i}./255;
    end
    figure
    b = barwitherr(Sems,Avgs);
%     title('ctrl internodes of stable designation from day 0-14')
    set(b(1),'FaceColor', colors{1});
    set(b(2),'FaceColor', colors{2});
    set(b(3),'FaceColor', colors{3});
    ylim([0 80])
    legend('0 neighbours','1 neighbour','2 neighbours', 'Location','northwest')
    xticklabels({'day 0' 'day 14'})
    ylabel('average myelin sheath length (\mum)')
    box off
    figQuality(gcf,gca,[4.5 4])
    
    stblB = [ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,1);...
        ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,2);...
        ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,3)];
    nbrs = [zeros(10,1);ones(10,1);ones(10,1).*2];
    [pStblNodeClassBslnCtrl,tbl,stats] = anova1(stblB,nbrs);
    pStblNodeClassBslnCtrl
%     multcompare(stats)
    
    stblL = [ALL.stats.NodeClassLnthAvg_last.ctrl(:,1);...
        ALL.stats.NodeClassLnthAvg_last.ctrl(:,2);...
        ALL.stats.NodeClassLnthAvg_last.ctrl(:,3)];
    nbrs = [zeros(10,1);ones(10,1);ones(10,1).*2];
    [pStblNodeClassLastCtrl,tbl,stats] = anova1(stblL,nbrs);
    pStblNodeClassLastCtrl
    
    [~,pStblBsln_Last] = ttest2(stblB,stblL)
%     multcompare(stats)
    %%
    Avgs = [mean(ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,4:6),'omitnan'); mean(ALL.stats.NodeClassLnthAvg_last.ctrl(:,4:6),'omitnan')]; 
    Sems = [calcSEM(ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,4:6)); calcSEM(ALL.stats.NodeClassLnthAvg_last.ctrl(:,4:6))]; 
    
    figure
    b = barwitherr(Sems,Avgs);
    hold on
%     title('ctrl dynamic internodes day 0-14')
    set(b(1),'FaceColor', colors{1});
    set(b(2),'FaceColor', colors{1});
    set(b(3),'FaceColor', colors{2});        
    legend('0->1 nbr','0->2 nbrs','1->2 nbrs','Location','northwest')
    xticklabels({'de novo' 'remyel'})
    ylabel('average myelin sheath length (\mum)')
    ylim([0 90])
    figQuality(gcf,gca,[4.5 4])
    hold off
    
    figure
    b = bar(Avgs);
    set(b(1),'FaceColor', colors{2});
    set(b(2),'FaceColor', colors{3});
    set(b(3),'FaceColor', colors{3}); 
    legend('0->1 nbr','0->2 nbrs','1->2 nbrs','Location','northwest')
    xticklabels({'de novo' 'remyel'})
    ylabel('average myelin sheath length (\mum)')
    ylim([0 90])
    figQuality(gcf,gca,[4.5 4])
    applyhatch_plusC(gcf,{makehatch_plus('\\4',5),makehatch_plus('\\4',5),makehatch_plus('\\4',5)},[colors{2};colors{3};colors{3}],[],500,2)
    
    dynamB = [ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,4);...
        ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,5);...
        ALL.stats.NodeClassLnthAvg_bsln.ctrl(:,6)];
    nbrs = [zeros(10,1);ones(10,1);ones(10,1).*2];
    [pStblNodeClassBslnCtrl,tbl,stats] = anova1(dynamB,nbrs);
    pStblNodeClassBslnCtrl
%     multcompare(stats)
    
    dynamL = [ALL.stats.NodeClassLnthAvg_last.ctrl(:,4);...
        ALL.stats.NodeClassLnthAvg_last.ctrl(:,5);...
        ALL.stats.NodeClassLnthAvg_last.ctrl(:,6)];
    nbrs = [zeros(10,1);ones(10,1);ones(10,1).*2];
    [pStblNodeClassLastCtrl,tbl,stats] = anova1(dynamL,nbrs);
    pStblNodeClassLastCtrl
    multcompare(stats)
    
    [~,pStblBsln_Last] = ttest2(dynamB,dynamL)
    
    %%
    % compare ctrl stable and dynamic ints at baseline
    [~,pBsln_StblDynam] = ttest2(stblB,dynamB)   
    % compare ctrl stable and dynamic ints at day 14
    [~,pLast_StblDynam] = ttest2(stblL,dynamL)
    
    avgCtrl = [mean(stblB,'omitnan') mean(dynamB,'omitnan') mean(stblL,'omitnan') mean(dynamL,'omitnan')];
    semCtrl = [calcSEM(stblB') calcSEM(dynamB') calcSEM(stblL') calcSEM(dynamL')];
    %%
    stblBcupr = [ALL.stats.NodeClassLnthAvg_bsln.cupr(:,1);...
            ALL.stats.NodeClassLnthAvg_bsln.cupr(:,2);...
            ALL.stats.NodeClassLnthAvg_bsln.cupr(:,3)];
    stblLcupr = [ALL.stats.NodeClassLnthAvg_last.cupr(:,1);...
            ALL.stats.NodeClassLnthAvg_last.cupr(:,2);...
            ALL.stats.NodeClassLnthAvg_last.cupr(:,3)];
    dynamBcupr = [ALL.stats.NodeClassLnthAvg_bsln.cupr(:,4);...
            ALL.stats.NodeClassLnthAvg_bsln.cupr(:,5);...
            ALL.stats.NodeClassLnthAvg_bsln.cupr(:,6)];
    dynamLcupr = [ALL.stats.NodeClassLnthAvg_last.cupr(:,4);...
            ALL.stats.NodeClassLnthAvg_last.cupr(:,5);...
            ALL.stats.NodeClassLnthAvg_last.cupr(:,6)];   
    avgCupr = [mean(stblBcupr,'omitnan') mean(dynamBcupr,'omitnan') mean(stblLcupr,'omitnan') mean(dynamLcupr,'omitnan')];
    semCupr = [calcSEM(stblBcupr') calcSEM(dynamBcupr') calcSEM(stblLcupr') calcSEM(dynamLcupr')];
    
    alldata = [stblB;stblL;dynamB;dynamL;stblBcupr;stblLcupr;dynamBcupr;dynamLcupr];
    types1={};types2={};types3={};types4={};types5={};types6={};types7={};types8={};
    types1(1:length(stblB)) = {'stblB'};
    types2(1:length(stblL)) = {'stblL'};
    types3(1:length(dynamB)) = {'dynamB'};    
    types4(1:length(dynamL)) = {'dynamL'};
    types5(1:length(stblBcupr)) = {'stblBcupr'};
    types6(1:length(stblLcupr)) = {'stblLcupr'};
    types7(1:length(dynamBcupr)) = {'dynamBcupr'};
    types8(1:length(dynamLcupr)) = {'dynamLcupr'};
    types = [types1,types2,types3,types4,types5,types6,types7,types8]';
    [p,tbl,stats] = anova1(alldata,types);
end