clc; clear

 pe = pyenv(Version='C:\Users\santo\AppData\Local\Programs\Python\Python310\python.exe');
% 
% %pysys = py.sys.path;
% 
% %pysys.append('C:\users\santo\appdata\local\packages\pythonsoftwarefoundation.python.3.11_qbz5n2kfra8p0\localcache\local-packages\python311\site-packages')
% %np = py.importlib.import_module('numpy');
% 
% 
% pyrunfile("main_EntropiaEDmatlab.py")
% 

pe = pyenv;
if pe.Status == "NotLoaded"
    [~,exepath] = system("where python");
    pe = pyenv('Version',exepath);
end