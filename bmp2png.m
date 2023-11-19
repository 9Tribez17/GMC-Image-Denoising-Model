Input_path = 'data2\';  
Output_path='result2\';
namelist = dir(strcat(Input_path,'*.bmp'));  %获得文件夹下所有的 .jpg图片
len = length(namelist);
for i = 1:len
    name=namelist(i).name;  %namelist(i).name; %这里获得的只是该路径下的文件名
    I=imread(strcat(Input_path, name));%图片完整的路径名

    imwrite(I,[Output_path,erase(name, ".bmp"),'.png']); %完整的图片存储的路径名  并将整形的数字 
end
