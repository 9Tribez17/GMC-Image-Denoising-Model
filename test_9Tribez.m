        clear all
        path ="data2/yacht.bmp";   
        Img = im2double(im2gray(imread(path)));
        
        hsize = 5;
        std = 7; % 5 %
        H = fspecial('gaussian',hsize,std);
        Bim = imfilter(Img,H,'circular','conv');
        %figure;imshow(Bim,[])
        sigma = 25/255;
        f = imnoise(Bim,'gaussian',0,sigma^2);
        
        Blu_g = f;
        
        %Ker_Size= [5 5];
        %Blu_Ker = fspecial('average', Ker_Size);

        %Sigma = 0.008;
        %Blu_g =  imfilter(Img, Blu_Ker,'circular','conv');%circular   replicate symmetric
        %  Lip_C = Est_Lip_Ker_Mat_C(ones(512),Blu_Ker,Ker_Size);
        psnr_bnim = psnr(Blu_g,Img);
        %Blu_g = Blu_g + Sigma*randn(size(Img));
        %figure(2); imshow(Blu_g,[]);
        %title(['PSNR = ' num2str(psnr_bnim)]);

        %I = Img;
        im = Img;
        Bnim = Blu_g;
        % n = size(I,1);
        miter = 500;
        %H = Blu_Ker; 
        Miter = 100;
        Bcho = 2;

        %disp('The ADMM for TGV regularization model is implemented');
        %disp('The optimal regularization parameters are nu1 = 0.01, nu0 = 0.03 for noiseL = 10');
        %disp('The optimal regularization parameters are nu1 = 0.01, nu0 = 0.01 for noiseL = 15');
        %disp('The optimal regularization parameters are nu1 = 0.012, nu0 =0.03 for Aerial');
        
        nu1list = [0.03];%[0.05:0.01:0.1];
                  %[0.02:0.01:0.06];
        nu2list = [0.03];
        result_TGV = [];
        [m1,n1] = size(H); 
        A = zeros(m1,n1);
        
        H1 = zeros(m1,n1);
        for i = 1 : m1
            A(i,:) = H(m1+1-i,:);
        end
        for j = 1 : n1
            H1(:, n1+1-j) = A(:, j);
        end
        
        for i = 1:length(nu1list)
            for j =1:length(nu2list)
                nu1 = nu1list(i);
                nu2 = nu2list(j);
                nu = [nu1,nu2];
                %nu = [0.01,0.03];
                %B = fft2(imfilter(Bnim,H1,nu,'circular','conv'));
                %B = fft2(imfilter(Bnim,H1,'circular','conv'));
                %if Bcho == 2
                [u,psnr_TGVL] = TGV2L2_ADMM(im,Bnim,H,nu,Miter);
%               else 
%               [u,psnr_TGVL] = TGV2L2_LADMM(im,Bnim,H,H1,alpha,delta,miter,Bcho);  
                psnr_TGV = psnr(u,im);
                result_TGV = [result_TGV;nu(1),nu(2),psnr_TGV]
                
            end
        end

 %%
 result_TV = [];
 lamlist   = [0.03];
 rholist   = [2];

%disp('The optimal regularization parameters are lam = 0.01, rho = 2 for noiseL = 10');
%disp('The optimal regularization parameters are lam = 0.01, rho = 5 for noiseL = 15');

 for i = 1:length(lamlist)
     for j =1:length(rholist)
         lam = lamlist(i);
         rho = rholist(j);
         [psnr_TV,imtv] = TV_ADMM(Img,f,H,lam,rho);
         result_TV = [result_TV;lam,rho,psnr_TV]
         
     end
 end

%%
result_GMC = [];
rho1list   = [2];
rho2list   = [60];
gammalist  = [0.03];
lambda1    = [0.03]; %正则化参数

%disp('The optimal regularization parameters are lam = 0.01, rho1 = 2,  rho2 = 30 for noiseL = 10');
%disp('The optimal regularization parameters are lam = 0.01, rho1 = 15, rho2 = 30 for noiseL = 15');
for i = 1:length(lambda1)
     for j = 1:length(rho1list)
        for m = 1:length(rho2list)
            for n = 1:length(gammalist)
                lam = lambda1(i);
                rho1 = rho1list(j);
                rho2 = rho2list(m);
                gamma = gammalist(n);
                [psnr_GMC, uad] = GMC_ADMM(Img,f,H,lam,rho1,rho2,gamma);
                result_GMC = [result_GMC;lam,rho1,rho2,gamma,psnr_GMC]
                
            end
        end
     end
end
Output_path='result2\';
imwrite(Img,[Output_path,'yacht_gray.png']);
imwrite(f,[Output_path,'yacht_noise.png']);
imwrite(u,[Output_path,'yacht_tgv.png']);
imwrite(imtv,[Output_path,'yacht_tv.png']);
imwrite(uad,[Output_path,'yacht_gmc.png']);


