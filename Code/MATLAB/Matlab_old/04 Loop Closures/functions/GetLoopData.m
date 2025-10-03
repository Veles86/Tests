function [loopData] = GetLoopData(loop,fullData)

n = length(loop);

loopData = [];

for i=1:n-1
    
    temp1 = fullData((fullData.A==loop(i)) & ((fullData.B==loop(i+1))),:);
    temp2 = fullData((fullData.A==loop(i+1)) & ((fullData.B==loop(i))),:);

    temp = [temp1;InvertRows(temp2)];

    loopData = [loopData;temp(1,:)];

end

    temp1 = fullData((fullData.A==loop(n)) & ((fullData.B==loop(1))),:);
    temp2 = fullData((fullData.A==loop(1)) & ((fullData.B==loop(n))),:);

    temp = [temp1;InvertRows(temp2)];

    loopData = [loopData;temp(1,:)];


end


function [invertedRows] = InvertRows(inputRows)


    invertedRows = inputRows;
    invertedRows(:,["A","B"]) = inputRows(:,["B","A"]);
    invertedRows.dH = -inputRows.dH;

    invertedRows(:,["codeA","typeA","nameA","groupA","rankA"]) = inputRows(:,["codeB","typeB","nameB","groupB","rankB"]);
    invertedRows(:,["codeB","typeB","nameB","groupB","rankB"]) = inputRows(:,["codeA","typeA","nameA","groupA","rankA"]);
    invertedRows(:,["Y1","X1","H1"]) = inputRows(:,["Y2","X2","H2"]);
    invertedRows(:,["Y2","X2","H2"]) = inputRows(:,["Y1","X1","H1"]);
    invertedRows{:,["dH_from_DB","dH_from_file"]} = inputRows{:,["dH_from_DB","dH_from_file"]};

end


