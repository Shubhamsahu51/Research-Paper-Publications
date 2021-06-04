clear all 
clc;
close all;
a1=imread('SAR_image_2.jpg');
a=rgb2gray(a1);
R=a;
a=imresize(a,[256 256]);
I=imtool(R)
whos;
[row col]=size(a);
a1=a;
figure(11) 
imshow(uint8(a));
 title('Original SAR img');

figure(12) 
imhist(a) ;
title('Original img Histogram');

a = imnoise(a,'speckle',0.02);

figure(1) 
imshow(uint8(a));
 title('Noisy SAR img');
figure(41) 
imhist(a) ;
title('Noisy img Histogram');
%  %Add multiplicative noise to the image
% a = imnoise(a,'speckle',0.01);
% figure(2),imshow(a);
% title('Noisy SAR img');

% DWT decomposition
[ NN,MM ]=size(a);
final=zeros(NN);
a1=double(a);
[ row col ]=size(a1);
[E1 F1 G1 H1]=dwt2(a1,'db1');
first=[E1 F1;G1 H1]
[I1 J1 K1 L1]=dwt2(E1,'db1');
second=[I1 J1 ;K1 L1];
final=[second F1;G1 H1];

figure(2)
imshow(uint8(first));
 title('first level DWT img');
figure(3)
imshow(uint8(final));
title('second level DWT img');
figure(4)
imshow(uint8(final));
title('second level DWT Coefficient img');


[row col]=size(I1);

%Apply Lee filter

OIm=myLee(I1);

%Apply Forust filter

[x y ]=size(I1);
I=I1;
K=1;
N=I;
for i=1:x
    for j=1:y                              
        if (i>1 & i<x & j>1 & j<y)
            mat(1)=I(i-1,j);
            mat(2)=I(i+1,j);
            mat(3)=I(i,j-1);
            mat(4)=I(i,j+1);
            d(1)=sqrt((i-(i-1))^2);
            d(2)=sqrt((i-(i+1))^2);
            d(3)=sqrt((j-(j-1))^2);
            d(4)=sqrt((j-(j+1))^2);
            mn=mean(mean(mat));
            c=mat-mn;
            c2=c.^2;
            c3=c/(c2+.0000001);
            Cs=0.25*sum(sum(c3));
            m(1)=exp(-K*Cs*d(1));
            m(2)=exp(-K*Cs*d(2));
            m(3)=exp(-K*Cs*d(3));
            m(4)=exp(-K*Cs*d(4));
            ms=sum(sum(m));
            mp=m/ms;
            N(i,j)=sum(sum(mp.*mat));                    
        end
     end
end
% ft=uint8(N);
ft = N;


% filtering with Kuan filter
[RI] = Kuanfilter(I1); 

% filtering with wiener filter
[WF] = wiener2(I1); 

% filtering with median filter
[MF] = medfilt2(I1,[3 3]); 

% filtering with meandi bivariant Shrikage Function usibg Forst
opt = 'gbl'; % Global threshold
T = 20;    % Threshold
[w1] = bishrink(I1,WF,T);

% finding inverse DWT
sx=size(final);

op2=idwt2(RI,J1,K1, L1 , 'db1', sx);
output2=idwt2(op2, F1,G1,H1, 'db1', sx);

op3=idwt2(ft,J1,K1, L1 , 'db1', sx);
output3=idwt2(op3, F1,G1,H1, 'db1', sx);

op4=idwt2(MF,J1,K1, L1 , 'db1', sx);
output4=idwt2(op4, F1,G1,H1, 'db1', sx);

op5=idwt2(WF,J1,K1, L1 , 'db1', sx);
output5=idwt2(op5, F1,G1,H1, 'db1', sx);

op6=idwt2(w1,J1,K1, L1 , 'db1', sx);
output6=idwt2(op6, F1,G1,H1, 'db1', sx);

op7=idwt2(OIm,J1,K1, L1 , 'db1', sx);
output7=idwt2(op7, F1,G1,H1, 'db1', sx);


figure(51)
imshow(uint8(output7));
title('inverse DWT Lee filtered image');
figure(52) 
imhist(uint8(output7) );
title('Histogram of Lee filtered image');

figure(5)
imshow(uint8(output2));
title('inverse DWT Kuan filtered image');
figure(45) 
imhist(uint8(output2)) ;
title('Histogram of Kuan filtered image');

figure(6)
imshow(uint8(output3));
title('inverse DWT fourst filtered image');
figure(46)
imhist(uint8(output3)) ;
title('Histogram of Frost filtered image');

 figure(7)
 imshow(uint8(output4));
title('inverse DWT Median filtered img');
figure(47)
imhist(uint8(output4)) ;
title('Histogram of median filtered image');


figure(8)
imshow(uint8(output5));
title('inverse DWT wiener filtered image');
figure(48)
imhist(uint8(output5)) ;
title('Histogram of wiener filtered image');


figure(9)
imshow(uint8(output6));
title('inverse DWT biv shk with wiener filter img');
figure(49)
imhist(uint8(output6)) ;
title('Histogram of biv shk filtered image');



% using wavelet fusion for better entropy
fuse_output = wfusimg(output6,output3,'db2',5,'mean','mean');

figure(10)
imshow(uint8(fuse_output));
title('wavelet fused ouptput img');
figure(410)
imhist(uint8(output6)) ;
title('Histogram of biv shk with fusion image');



 t1=0; t2=0; t3=0; t4=0; t5=0; t6=0; t7=0; t8=0;

for x=1:1:row
    for y=1:1:col
        Erms_8=(a1(x,y)-output7(x,y))^2;    %Kaun filtered MSE
        E8=t8+ Erms_8;
        t8=E8;        
               
        Erms_1=(a1(x,y)-output2(x,y))^2;    %Kaun filtered MSE
        E1=t1+ Erms_1;
        t1=E1;        
        
        Erms_2=(a1(x,y)-output3(x,y))^2;    %Fourst filtered MSE
        E2=t2+ Erms_2;
        t2=E2;        
        
        Erms_3=(a1(x,y)-output4(x,y))^2;    %Median filtered MSE
        E3=t3+ Erms_3;
        t3=E3;        
        
        Erms_4=(a1(x,y)-output5(x,y))^2;    %Wiener filtered MSE
        E4=t4+ Erms_4;
        t4=E4;        
        
        Erms_5=(a1(x,y)-output6(x,y))^2;    %Bivt_shk_forst filtered MSE
        E5=t5+ Erms_5;
        t5=E5;        
        
        Erms_7=(a1(x,y)-fuse_output(x,y))^2;    %Bivt_shk_forst filtered MSE
        E7=t7+ Erms_7;
        t7=E7;        
        
        or=(a1(x,y))^2;
        E6=t6+ or;
        t6=E6;
              
    end
end
C=row*col;
      Er_8=sqrt(t8/C);    %MSE Between Original and KUan image
     snr_8=sqrt(t6/t8);   %SNR Between Original and Kuan image
  
     Er_1=sqrt(t1/C);    %MSE Between Original and KUan image
     snr_1=sqrt(t6/t1);   %SNR Between Original and Kuan image
  
       Er_2=sqrt(t2/C);    %MSE Between Original and Forust image
     snr_2=sqrt(t6/t2);   %SNR Between Original and Forust image
  
     Er_3=sqrt(t3/C);    %MSE Between Original and Median image
     snr_3=sqrt(t6/t3);   %SNR Between Original and Median image
  
     Er_4=sqrt(t4/C);    %MSE Between Original and Wiener image
     snr_4=sqrt(t6/t4);   %SNR Between Original and Wiener image
  
     Er_5=sqrt(t5/C);    %MSE Between Original and Bvt shk image
     snr_5=sqrt(t6/t5);   %SNR Between Original and bvt shk image
     
    Er_7=sqrt(t7/C);    %MSE Between Original and Bvt shk image
     snr_7=sqrt(t6/t7);   %SNR Between Original and bvt shk image
     
     
     Mean_Breightness=[ mean(mean(uint8(a1))) mean(mean(uint8(output7)))  mean(mean(uint8(output2)))  mean(mean(uint8(output3)))  mean(mean(uint8(output4)))  mean(mean(uint8(output5)))  mean(mean(uint8(output6)))  mean(mean(uint8(fuse_output)))];          
     
     MSE=[ Er_8 Er_1 Er_2 Er_3 Er_4 Er_5 Er_7];
     PSNR=[10*log(snr_8) 10*log(snr_1) 10*log(snr_2) 10*log(snr_3) 10*log(snr_4) 10*log(snr_5) 10*log(snr_7)];
     
     Entropy=[entropy(uint8(a1)) entropy(uint8(output7)) entropy(uint8(output2)) entropy(uint8(output3)) entropy(uint8(output4)) entropy(uint8(output5)) entropy(uint8(output6)) entropy(uint8(fuse_output))];
  
     bvt_anls=[Er_5  snr_5  entropy(w1) entropy(fuse_output)]; 
     
     
     