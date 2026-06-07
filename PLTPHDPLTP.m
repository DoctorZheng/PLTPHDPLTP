function  [PLTP_value,HDPLTP_value] = PLTPHDPLTP(I)
[m,n,h] = size(I);
if h == 3
    I = rgb2gray(I);
end
I = double(I);

thresh1 = 4;
thresh0 = 1;
thresh2 = 1;
v = mean2(I);

pow3  = 3.^(9:-1:0);
pow3d = 3.^(7:-1:0);

uMat = conv2(I, ones(5)/25, 'same');

%% ===================== PLTP =====================
H = m-4;
W = n-4;
PLTP_value = zeros(H, W, 'double');
PLTP_data  = zeros(H*W, 10, 'double');  
idx = 1;

for i = 3:m-2
    i2s = i-2:i+2;
    for j = 3:n-2
        u = uMat(i,j);
        
        Ii = I(i,j);
        out_novertex = (I(i+2,j-2)+I(i+2,j)+I(i+2,j+2)+I(i+1,j-2)+I(i+1,j+2)...
                      + I(i-1,j-2)+I(i-1,j+2)+I(i-2,j-2)+I(i-2,j)+I(i-2,j+2)) * 0.1;
                  
        in_novertex  = (I(i-1,j-1)+I(i-1,j+1)+I(i,j-1)+I(i,j+1)+I(i+1,j)) * 0.2;
        
        out_vertex   = (I(i-2,j-1)+I(i-2,j+1)+I(i,j+2)+I(i+2,j+1)+I(i+2,j-1)+I(i,j-2))/6;
        in_vertex    = (I(i-1,j)+I(i+1,j+1)+I(i+1,j-1))/3;

        vec0 = [I(i-2,j-1) I(i-2,j+1) I(i,j+2) I(i+2,j+1) I(i+2,j-1) I(i,j-2) out_vertex out_novertex];
        vec1 = [I(i-1,j) I(i-1,j) I(i+1,j+1) I(i+1,j+1) I(i+1,j-1) I(i+1,j-1) in_vertex in_novertex];
        diff0 = vec0 - vec1;
        diff1 = vec1 - Ii;
        
        neighbor01 = diff0 > thresh1;
        neighbor02 = diff0 < -thresh1;
        neighbor11 = diff1 > thresh1;
        neighbor12 = diff1 < -thresh1;

        neighbor0 = neighbor01 - neighbor02;
        neighbor1 = neighbor11 - neighbor12;
        temp = neighbor0 + neighbor1;
        neighbor = (temp > 0) + (temp >= 0);

        n00 = [(Ii-v>thresh0), (Ii-u>thresh0)];
        n01 = [(Ii-v<-thresh0),(Ii-u<-thresh0)];
        n0 = n00 - n01 + 1;

  
        data = [n0, neighbor];
        PLTP_data(idx,:) = data;
        idx = idx + 1;

        PLTP_value(i-2,j-2) = (data * pow3') / 230;
    end
end

%% ===================== HDPLTP  =====================
[m1,n1] = size(PLTP_value);
Hd = m1-10;
Wd = n1-10;
HDPLTP_value = zeros(Hd, Wd, 'double');

pos = [-5 -5; -5 0; -5 5; 0 5; 5 5; 5 0; 5 -5; 0 -5];

for i = 6:m1-5
    ji = (i-1)*n1;
    for j = 6:n1-5
        c = PLTP_data(ji + j, :);
        td = zeros(1,8);
        
        for k = 1:8
            ni = i + pos(k,1);
            nj = j + pos(k,2);
            td(k) = sum(c ~= PLTP_data((ni-1)*n1 + nj, :));
        end
        
        tdmean = mean(td);
        tdd = td > tdmean + thresh2;
        tdu = td >= tdmean - thresh2;
        code = tdd + tdu;
        
        HDPLTP_value(i-5,j-5) = sum(code.*pow3d) / 25;
    end
end

PLTP_value   = int32(PLTP_value);
HDPLTP_value = int32(HDPLTP_value);
