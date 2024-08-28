close all;
clear all;

file_path =  '.\fig-test\1\';% Image folder path，前5个假，后5个真
cc=colormap(lines(100));
img_path_list = dir(strcat(file_path,'*.bmp'));
img_num = length(img_path_list);
NUMD0=zeros(img_num,10);
NUMD=zeros(img_num,10);
for ii = 1:10     
    I=imread(strcat('fig-test\',strcat(num2str(ii),'.bmp')));
    image_size=size(I);
    dimension=numel(image_size);
    if dimension~=2
       IT1 = rgb2gray(I);
    else
       IT1 = (I);             
    end
    [H W]=size(IT1);  
    
    [PLTPv1,HDPLTPv1]=PLTPHDPLTP(IT1);  

    PLTPv1=PLTPv1(1+5:H-4-5,1+5:W-4-5);

    %test
    step=4;
    file_path =  strcat(strcat('.\',num2str(ii)),'\');
    if img_num > 0 
        Dis=zeros(floor(H*W/(step*step)),img_num);
        for jj = 1:img_num 
            image_name = img_path_list(jj).name;
            I1 =  imread(strcat(file_path,image_name));
            
            image_size=size(I1);
            dimension=numel(image_size);
            if dimension~=2
               IT2 = rgb2gray(I1);
            else
               IT2 = (I1);                
            end                        
            [PLTPv2,HDPLTPv2]=PLTPHDPLTP(IT2);    
            PLTPv2=PLTPv2(1+5:H-4-5,1+5:W-4-5);
            
            Dnum=0;  
            nuuu=1;
            sized=5;
            for i = sized+1 : step : H-14-sized
                for j = sized+1 : step : W-14-sized                                         
                    PLTPv1IT=PLTPv1(i-sized:i+sized, j-sized:j+sized);
                    HDPLTPv1IT=HDPLTPv1(i-sized:i+sized, j-sized:j+sized);
                    PLTPv2IT=PLTPv2(i-sized:i+sized, j-sized:j+sized);
                    HDPLTPv2IT=HDPLTPv2(i-sized:i+sized, j-sized:j+sized);
                    bins= 256; 
                    result1(1,:)=hist(PLTPv1IT(:),0:(bins-1));
                    result1(2,:)=hist(HDPLTPv1IT(:),0:(bins-1));
                    result1(1,:)=result1(1,:)/sum(result1(1,:));
                    result1(2,:)=result1(2,:)/sum(result1(2,:));
                    result2(1,:)=hist(PLTPv2IT(:),0:(bins-1));
                    result2(2,:)=hist(HDPLTPv2IT(:),0:(bins-1));
                    result2(1,:)=result2(1,:)/sum(result2(1,:));
                    result2(2,:)=result2(2,:)/sum(result2(2,:));
                    Dis(nuuu,jj)=pdist2(result1(:)',result2(:)','euclidean');
                    if  Dis(nuuu,jj)>0.19 
                        Dnum=Dnum+1;
                    end 
                    nuuu=nuuu+1;
                end
            end                   
            NUMD0(jj,ii)=Dnum;
        end
    end  
end


NUMD1=NUMD0/nuuu;
TT=0.84;   
AA=NUMD1;
AA(NUMD1<TT)=1;%true;   AA(6:10,:) is positive
AA(NUMD1>=TT)=0;%false;   AA(1:5,:) is negative
AA


