function handleNewFig = copyFig(handleOriginalFig)


if nargin == 0
    handleOriginalFig=gcf;
end
handleNewFig=figure;
objects=allchild(handleOriginalFig);
copyobj(get(handleOriginalFig,'children'),handleNewFig);