function [coordinatesTable] = GetCoordinatesFromMergedTable(InputMergedTable)

n = height(InputMergedTable);

coordinatesTable = table('Size',[2*n 8],'VariableTypes',{'double','string','string','double','double','double','double','double'},'VariableNames',["code","type","name","group","rank","Y","X","H"]);
coordinatesTable(1:n,:) = [InputMergedTable(:,15:19),InputMergedTable(:,25:27)];
coordinatesTable(n+1:end,:) = [InputMergedTable(:,20:24),InputMergedTable(:,28:30)];

rows = isnan(coordinatesTable.code);
coordinatesTable(rows,:) = [];

rows = isnan(coordinatesTable.rank);
coordinatesTable.rank(rows) = 9999;
rows = isnan(coordinatesTable.group);
coordinatesTable.group(rows) = 9999;
rows = isnan(coordinatesTable.H);
coordinatesTable.H(rows) = 9999;
coordinatesTable = unique(coordinatesTable,'rows');


rows = coordinatesTable.rank==9999;
coordinatesTable.rank(rows) = NaN;
rows = coordinatesTable.group==9999;
coordinatesTable.group(rows) = NaN;
rows = coordinatesTable.H==9999;
coordinatesTable.H(rows) = NaN;

coordinatesTable = sortrows(coordinatesTable,1);

end



