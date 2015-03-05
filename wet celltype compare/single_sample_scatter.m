function single_sample_scatter(measured_profile, predicted_profile)
    transparentScatter(measured_profile, predicted_profile,0.005,0.2);
%     scatter(measured_profile, predicted_profile,'.')
    title(sprintf('spearman - %g', corr(measured_profile, predicted_profile,'type','spearman')));
    xlabel('measured profile');
    ylabel('predicted profile');
    xLimits = get(gca,'XLim');xLimits(1) = 0;
    yLimits = get(gca,'YLim');yLimits(1) = 0;
    xlim(xLimits); ylim(yLimits);
%     axis equal;
    max_x = ceil(max(measured_profile));
    max_y = ceil(max(predicted_profile));
    set(gca,'ytick',0:max_y)
    set(gca,'xtick',0:max_x)
end