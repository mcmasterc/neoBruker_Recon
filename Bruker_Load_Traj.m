function traj = Bruker_Load_Traj(file_name);

%Get the file
if nargin>0
      filename= file_name;  
else
    [filename,path]=uigetfile('*','Select file');
    cd(path)
end

% open data
fd=fopen(char('traj'),'r','l'); %n, b, l, s, a
fseek(fd,0,'bof');  
raw=fread(fd,inf,'double');
fclose(fd);
traj = reshape(raw,3,[]);

end
