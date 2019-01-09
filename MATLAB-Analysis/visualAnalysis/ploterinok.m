figure
hold on
for ii = 1:length(s.waves_03)-1
    plot(s.waves_03(ii,17:47),'color',[0.7 0.7 0.7]);
end
plot(mean(s.waves_03(:,17:47),1),'color',[0 0 0]);