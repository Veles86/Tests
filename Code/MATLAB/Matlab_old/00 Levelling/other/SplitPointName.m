function [ number, letter ] = SplitPointName( fullPointName )

if ~isa(fullPointName,'char')
    fullPointName = char(fullPointName);
end

ascii = double(fullPointName);
digits = (ascii > 47 & ascii < 58);

if (sum(digits)==length(digits))
    if length(digits)<=4
        number = string(fullPointName);
        letter = string('');
    else
        number = string(fullPointName(1:4));
        letter = string(fullPointName(5:end));
    end
elseif sum(digits)==0
    number = string('');
    letter = string(fullPointName);
else
    index = find(digits==0,1);
    if index<=5
        number = string((fullPointName(1:index-1)));
        letter = string(fullPointName(index:end));
    else
        number = string((fullPointName(1:4)));
        letter = string(fullPointName(5:end));
    end
end

end