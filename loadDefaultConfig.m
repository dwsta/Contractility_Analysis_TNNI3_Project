function [cfg_data] = loadDefaultConfig()
configfile = 'default_config.json';
fid = fopen(configfile);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
cfg_data = jsondecode(str);
end