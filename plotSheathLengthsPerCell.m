function plotSheathLengthsPerCell(data,ALL)
cond = {'ctrl','bsln','L2_3'};
figure
hold on
x = 0;
for c = 1:length(cond)
   cells = fieldnames(ALL.(cond{c}));
    for k = 1:length(cells)
        days = fieldnames(ALL.(cond{c}).(cells{k}).raw);
        if contains(days{end},{'d12','d14','d15'}) && isfield(data.(cond{c}),cells{k})
            lnths = ALL.(cond{c}).(cells{k}).raw.(days{end}).Sheath_Info.intLengths;
            brgstatus = ALL.(cond{c}).(cells{k}).raw.(days{end}).Sheath_Info.isbridgecnctd;
            brgstatus(lnths==0) = [];
            lnths(lnths==0) = [];            
            x = x+1;
%             plot(x*ones(1,length(lnths(brgstatus==0)))',lnths(brgstatus==0),'ko');
%             plot(x*ones(1,length(data.(cond{c}).(cells{k}).lengthSums)),data.(cond{c}).(cells{k}).lengthSums,'ro','MarkerFaceColor','r')
            plotSpread(lnths(brgstatus==0),'xValues',x,'distributionColors','k');
            if any(data.(cond{c}).(cells{k}).lengthSums)
                y = data.(cond{c}).(cells{k}).lengthSums;
                y = y(~isnan(y));
                plot(x*ones(1,length(y)),y,'ro','MarkerSize',4)
            end
        end
    end
end