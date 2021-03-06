function _ = write_text_file(data, title_)
    % Write a text file from the matrix "data"
    file = fopen(strcat("./autogenerated_data/", title_, ".txt"), 'w');

    % 9 strings with one space between each on the first line
    fprintf(file, '%1s %1s %1s %1s %1s %1s %1s %1s %1s\r\n', "t", "Wilson", "Newmark1", "Newmark2", "Newmark3",
    "Central Difference", "Constant Approximation", "Linear Approximation", "Analytical");

    % 9 numbers with one space between each for every row in the matrix "data"
    fprintf(file, '%1f %1f %1f %1f %1f %1f %1f %1f %1f\r\n', data);
    fclose(file);
end
