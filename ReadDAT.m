function B = ReadDAT(image_size,data_path)
%% function B = ReadDAT(image_size, data_path)
%read the '.dat' file
%%input :
%    image_size  : image size
%    data_path   : '.dat' file
%%output:
%    B           : the output label image where the value at each pixel stands for its superpixel label 

row = image_size(1);
colomn = image_size(2);
fid = fopen(data_path,'r');
A = fread(fid, row * colomn, 'uint32')';
A = A + 1;
B = reshape(A,[colomn, row]);
B = B';
fclose(fid);