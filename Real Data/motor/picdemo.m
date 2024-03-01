h=gca;
YourYticklabel=cell(size(h.YDisplayLabels));
[YourYticklabel{:}]=deal('');
h.YDisplayLabels=YourYticklabel;
YourYticklabel=cell(size(h.XDisplayLabels));
[YourYticklabel{:}]=deal('');
h.XDisplayLabels=YourYticklabel;