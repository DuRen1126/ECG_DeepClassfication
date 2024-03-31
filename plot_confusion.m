function plot_confusion(confusion)
    figure();
    heatmap({'L', 'N', 'R', 'V'}, {'L', 'N', 'R', 'V'},confusion);
    xlabel('Predicted labels');
    ylabel('True labels');
    title('Confusion Matrix');
end
