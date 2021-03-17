function [m,e,h1,h2] = make_my_bar_nan(position, data, face_color, line_color, bar_width, face_alpha)

if nargin<6
    face_alpha=1;
end

m = nanmean(data);
e = nanstd(data)/sqrt(   nnz(not(isnan(data)))      );
h1=bar(position, m,'FaceColor', face_color, 'EdgeColor', line_color, 'LineWidth', .6, 'BarWidth',bar_width, 'FaceAlpha', face_alpha);
h2=plot([position; position], m+[-e; +e], 'Color', line_color, 'LineWidth', 1.5);

end
