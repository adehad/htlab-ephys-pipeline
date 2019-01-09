saveFig = 0;

d = n+n1;
for ii = 1:length(b)
    f{ii} = cat(1,b{ii},b1{ii});
    h{ii} = cat(1,q{ii},q1{ii});
end

x = '02';
e = linspace(0,1,50);
c = e(1:end-1)+e(2)/2;
thresh = 0.1;
savefig = 1;

u = 0;
% for kk = 0:2:6
%     figure
%     hold on
%     for jj = [1 9 17 25 31 39]+kk
%         u = u + 1;
%         plot(c,d(u,:)/(48+49))
%     end
%     set(gcf,'color','w');
%     set(gcf,'renderer','Painters');
%     xlabel('Normalised distance along trajectory line'); ylabel('Spike count per trial');
%     legend()
%     hold off
%     ylim([0,1])
%     
%     % save figures
%     if saveFig == 2
%         export_fig(sprintf('%s_trajHist_unit_%s_post.eps',opt.preName,s.clusters(ii)))
%         saveas(gcf, join([opt.preName '_trajHist_unit_' x '_' num2str(kk/2)],''), 'fig');
%     elseif saveFig
%         saveas(gcf, join([opt.preName '_trajHist_unit_' x '_' num2str(kk/2)],''), 'epsc');
%         saveas(gcf, join([opt.preName '_trajHist_unit_' x '_' num2str(kk/2)],''), 'fig');
%     end
% end

for ii = 1:length(b)
    tBins = find(d(ii,:)/97 > thresh);
    tSpikes = ismember(f{ii},tBins);
    tLoc = h{ii}(tSpikes);
    figure
    scatter(tLoc,zeros(size(tLoc)),2*ones(size(tLoc)),'r');
    hold on
    plot([0 1],[0 0],'color',[0.5 0.5 0.5])
    xlim([0 1])
    saveas(gcf, join([opt.preName '_trajHist_unit_' x '_line_' num2str(ii)],''), 'epsc')
end
