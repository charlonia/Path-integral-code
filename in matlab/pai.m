clear
clc
tic
%% 上色
a = imread('rgb.png');
mapea=ones(256,3);
for k =1:3
for j=1:256
    mapea(257-j,k)=double(a(round(j/256*490)+20,1,k))/256;
end
end
%%
grid = 1 ;%格点数目
grid_div = 0.01;% 格点精准度
alen = grid/grid_div;
time_div = 0.002;%时间间隔
time_all = 0.006;%总的演化时间
d = 0.4     %狭缝宽度  
v=10; % 设置电子入射速度
hb=1;
m=1;
A = zeros(grid/grid_div,grid/grid_div,2);%储存在任意时间任意地点电子状态的cell
mkdir(string(grid_div)+'r_con'+string(time_div)+'d='+string(d))
mkdir('b '+string(grid_div)+'r_con'+string(time_div)+'d='+string(d))

%%
%单缝衍射
Numx=meshgrid(1:1:grid/grid_div);
Numy=meshgrid(1:1:grid/grid_div)';
z = -i*grid_div^2/(2*time_div)*m/hb;

%/()
% 输入平面波 本身方程为exp(i*p*x/h_bar)，这里设置m和h_bar 为1
f = @(x) exp(i*m*v*x/hb);
for j=1:grid/2/grid_div-1
    A(:,j,1)=A(:,j,1)+f(j*grid_div);
end
% C = A(:,:,1).*conj(A(:,:,1));
% A(:,:,1)=A(:,:,1)/sum(sum(C));
L = real(A(:,:,1));
L = L-min(min(L));
image(256*L/(max(max(L))));



%% 开始计算迭代A(:,:,1)记录上一次的结果，A(:,:,2)记录由A(:,:,1)得到的结果，再赋值回去
for j =2:time_all/time_div
%     for l=1:grid/2/grid_div-1
%     A(:,l,1)=A(:,l,1)+f(l*grid_div-v*time_div*j);
%     end
%     result = sum(sum(A(:,:,1)));
%     result = (conj(result)*result)^0.5
%     director(j) = result;
     %% 传播子计算
     %左半边
    for kx=1:grid/grid_div/2-1
     for ky=1:grid/grid_div
        G=exp(z.*((Numx-kx).^2+(Numy-ky).^2)).*A(:,:,1).*grid_div.^2;
        G = sum(sum(G(1:grid/grid_div,1:grid/grid_div/2-1)));
        a=G;
        for lx=grid/grid_div/2:grid/grid_div-1
            for ly=max(round(((grid-d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),1):min(round(((grid+d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),alen)
                a=a+exp(z*((lx-kx)^2+(ky-ly)^2))*A(ly,lx,1)*grid_div.^2;
            end
        end
        A(ky,kx,2)=a;
     end
    end
    % 右半边
    for kx=grid/grid_div/2+1:grid/grid_div
    for ky=1:grid/grid_div
        a =0;
        G=exp(z.*((Numx-kx).^2+(Numy-ky).^2)).*A(:,:,1).*grid_div.^2;
        G = sum(sum(G(1:grid/grid_div,grid/grid_div/2+1:grid/grid_div-1)));
        a=G;
        for lx=1:grid/grid_div/2
            for ly=max(round(((grid-d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),1):min(round(((grid+d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),alen)
                a=a+exp(z*((lx-kx)^2+(ky-ly)^2))*A(ly,lx,1)*grid_div.^2;
            end
        end
        if kx==grid/grid_div
          a=a+A(ky,kx,1);
        end
        A(ky,kx,2)=a;
      % 屏幕吸收光子并且不再传播，屏幕位置就在grid/grid_div处
    end
    end
    kx=grid/grid_div/2;
    %中间 
    for ky=1:grid/grid_div
        G=exp(z.*((Numx-kx).^2+(Numy-ky).^2)).*A(:,:,1).*grid_div.^2;
        G = sum(sum(G(1:grid/grid_div,1:grid/grid_div-1)));
        a=G;
        A(ky,kx,2)=a;
    end
A(1:round((grid-d)/grid_div/2),grid/grid_div/2,1)=0;
A(round((grid+d)/grid_div/2+1):grid/grid_div,grid/grid_div/2,1)=0;
A(:,:,1)=A(:,:,2); 
% 强制归一化
C = A(:,:,1).*conj(A(:,:,1));
% A(:,1:alen-1,1)=A(:,1:alen-1,1)/sum(sum(C(:,1:alen-1)));
%%
%可视化
L = real(A(:,:,1));
L = L-min(min(L));
image(256*L/(max(max(L))));
colorbar
imwrite(round(256*L/(max(max(L)))),mapea,string(grid_div)+'r_con'+string(time_div)+'d='+string(d)+'\phase_'+string(j)+'.tif');
imwrite(round(256*L/(max(max(L)))),'b '+string(grid_div)+'r_con'+string(time_div)+'d='+string(d)+'\phase_'+string(j)+'.tif');
image(256*C/(max(max(C))));
imwrite(round(256*C/(max(max(C)))),mapea,string(grid_div)+'r_con'+string(time_div)+'d='+string(d)+'\possi_'+string(j)+'.tif');
imwrite(round(256*C/(max(max(C)))),'b '+string(grid_div)+'r_con'+string(time_div)+'d='+string(d)+'\possi_'+string(j)+'.tif');
csvwrite(string(grid_div)+'r_con'+string(time_div)+'d='+string(d)+'\drector.csv',director)
j
end
toc
%% 理论曲线
% hef = @(x) (sin(m.*v.*d.*x.*grid_div./hb)./(m.*v.*d.*x.*grid_div./hb)).^2;
% Z = [alen/2:-1:1,1:1:alen/2];
% Y = hef(Z);
% plot(1:1:alen,hef(Z));
%% 狭缝测试

% kx =51 ;
% ky =50;
% 
%   for ly=1:grid/grid_div
%          for lx=grid/grid_div/2+1:grid/grid_div
%              A(ly,lx,1)=256;
%          end
%         end
%         for lx=1:grid/grid_div/2
%             for ly=max(round(((grid-d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),1):min(round(((grid+d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),alen)
%                A(ly,lx,1)=256;
%             end
%         end
        
% for ly=1:grid/grid_div
%          for lx=1:grid/grid_div/2-1
%              A(ly,lx,1)=256;
%          end
%         end
%         for lx=grid/grid_div/2:grid/grid_div-1
%             for ly=max(round(((grid-d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),1):min(round(((grid+d)/grid_div/2-ky)/(grid/grid_div/2-kx)*(lx-kx)+ky),alen)
%                A(ly,lx,1)=256;
%             end
%         end
%  A(ky,kx,1)=0;
%  image(A(:,:,1))



%         for lx=1:grid/grid_div/2-1
%             for ly=1:grid/grid_div
%              a=a+exp(z*((lx-kx)^2+(ky-ly)^2))*A(ly,lx,1)*grid_div;
%           end
%         end