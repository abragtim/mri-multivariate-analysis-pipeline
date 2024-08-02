function [json] = readjson(path)
    file = fileread(path);
    json = jsondecode(file);
end
