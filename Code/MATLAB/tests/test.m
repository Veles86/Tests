% Initialize an empty string
concatenatedString = "";

% Loop to add codes
for i = 1:5
    % Create a code string for the current iteration
    codeString = sprintf("%d", i);

    % Concatenate the new code string to the existing string
    concatenatedString = concatenatedString + codeString + ","; % Adding a space for separation
end

%Remove the trailing comma and space if it exists
if endsWith(concatenatedString, ',')
    concatenatedString = extractBefore(concatenatedString, length(char(concatenatedString)));
end

% Display the concatenated string
disp('Concatenated String:');
disp(concatenatedString);