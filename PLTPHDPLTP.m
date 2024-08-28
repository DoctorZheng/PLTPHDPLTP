function  [PLTP_value,HDPLTP_value] = PLTPHDPLTP(I)
[m,n,h] = size(I);
if h==3
    I =  rgb2gray(I);
end
I=double(I);

PLTP_value = double(zeros([m-4 n-4]));
HDPLTP_value = double(zeros([m-10-4 n-10-4]));
PLTP = cell(m-4,n-4);
HDPLTP = cell(m-10-4,n-10-4);
thresh1=4;   
thresh2 = 1;
v = mean2(I);
for i = 3:m-2
    for j = 3:n-2
        u = mean2(I(i-2:i+2,j-2:j+2)); 
        out_novertex = (I(i+2,j-2)+I(i+2,j)+I(i+2,j+2)+I(i+1,j-2)+I(i+1,j+2)+I(i-1,j-2)+I(i-1,j+2)+I(i-2,j-2)+I(i-2,j)+I(i-2,j+2))/10.0;
        in_novertex = (I(i-1,j-1)+I(i-1,j+1)+I(i,j-1)+I(i,j+1)+I(i+1,j))/5.0;
        out_vertex = (I(i-2,j-1)+I(i-2,j+1)+I(i,j+2)+I(i+2,j+1)+I(i+2,j-1)+I(i,j-2))/6.0;
        in_vertex = (I(i-1,j)+I(i+1,j+1)+I(i+1,j-1))/3.0;  
        neighbor01 = ([I(i-2,j-1) I(i-2,j+1) I(i,j+2) I(i+2,j+1) I(i+2,j-1) I(i,j-2) out_vertex out_novertex] - [I(i-1,j) I(i-1,j) I(i+1,j+1) I(i+1,j+1) I(i+1,j-1) I(i+1,j-1) in_vertex in_novertex])>thresh1;
        neighbor02 = ([I(i-2,j-1) I(i-2,j+1) I(i,j+2) I(i+2,j+1) I(i+2,j-1) I(i,j-2) out_vertex out_novertex] - [I(i-1,j) I(i-1,j) I(i+1,j+1) I(i+1,j+1) I(i+1,j-1) I(i+1,j-1) in_vertex in_novertex])<-thresh1;
        neighbor11 = ([I(i-1,j) I(i-1,j) I(i+1,j+1) I(i+1,j+1) I(i+1,j-1) I(i+1,j-1) in_vertex in_novertex] - I(i,j))>thresh1;
        neighbor12 = ([I(i-1,j) I(i-1,j) I(i+1,j+1) I(i+1,j+1) I(i+1,j-1) I(i+1,j-1) in_vertex in_novertex] - I(i,j))<-thresh1;        
        neighbor0 = neighbor01 - neighbor02;
        neighbor1 = neighbor11 - neighbor12;        
        neighbor_u = (neighbor0 + neighbor1) > 0;
        neighbor_d = (neighbor0 + neighbor1) >= 0;
        neighbor= neighbor_u + neighbor_d;
        
        n00=  [(I(i,j)- v > thresh2) (I(i,j)- u > thresh2)];
        n01=  [(I(i,j)- v < -thresh2) (I(i,j)- u < -thresh2)];
        n0 = n00 - n01 + 1;
        PLTP(i-2,j-2)={[n0 neighbor]};
        pixel = 0;
        for k = 1:8
            pixel = pixel + neighbor(1,k) * 3^(8-k);
        end
        pixel = pixel + n0(1) * 3^(8)+ n0(2) * 3^(9);
        PLTP_value(i-2,j-2) = double(pixel)/230;
    end
end

m1=m-4;
n1=n-4;
thresh3 = 1;
td = double(zeros([1 8]));
for i = 6:m1-5
    for j = 6:n1-5
        neighbord = [PLTP(i-5,j-5) PLTP(i-5,j) PLTP(i-5,j+5) PLTP(i,j+5) PLTP(i+5,j+5) PLTP(i+5,j) PLTP(i+5,j-5) PLTP(i,j-5)];
        for k = 1:8
            td(k) = sum((cell2mat(PLTP(i,j)) - cell2mat(neighbord(k)))~=0);%pdist([cell2mat(PLTP(i,j));cell2mat(neighbord(k))], 'hamming') * 10;
        end
        tdmean=mean(td);
        tdd=(td-tdmean-thresh3)>0;
        tdu=(td-tdmean+thresh3)>=0;
        HDPLTP(i-5,j-5)={tdd+tdu};
        pixel_d = 0;
        for k = 1:8
            pixel_d = pixel_d + (tdd(k)+tdu(k)) * 3^(8-k);
        end        
        HDPLTP_value(i-5,j-5) = double(pixel_d)/25;
    end
end
PLTP_value=int32(PLTP_value); 
HDPLTP_value=int32(HDPLTP_value); 
% bins= 256;%8*(8-1) + 3;  
% result(1,:)=hist(PLTP_value(:),0:(bins-1));%/ numel(lbpIu);
% result(2,:)=hist(HDPLTP_value(:),0:(bins-1));%/ numel(lbpIu);
% result(1,:)=result(1,:)/sum(result(1,:));
% result(2,:)=result(2,:)/sum(result(2,:));