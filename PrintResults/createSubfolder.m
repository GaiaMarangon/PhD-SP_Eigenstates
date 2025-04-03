function createSubfolder(folderName)
    % checks if the folder folderName exists, otherwise creates it
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
end