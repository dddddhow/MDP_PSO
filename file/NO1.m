clc,clear;
%叠后地震数据（sei）读取及参数设置
nt=8000;
nx=200;
dt=0.0004;
sei=zeros(nt,nx);
fp=fopen('record','rb');
sei=fread(fp,[nt nx],'float32');
fclose(fp);
% imagesc(sei)

% %一、计三瞬剖面(振幅剖面amp、相位剖面ang、频率剖面fre)
% hsei=zeros(nt,nx);
% hsei=hilbert(sei);
% amp=sqrt(real(hsei).*real(hsei)+imag(hsei).*imag(hsei));
% figure(1)
% imagesc(amp),xlabel('道数/道'),ylabel('时间'),title('振幅剖面');
% ang=atan(imag(hsei)./real(hsei));
% figure(2)
% imagesc(ang),xlabel('道数/道'),ylabel('时间'),title('相位剖面');
% fre=diff(ang,1)/(2*pi);
% figure(3)
% imagesc(fre),xlabel('道数/道'),ylabel('时间'),title('频率剖面');


hsei = zeros(nt,nx);

for i=1:1:nt
    if(sei(i,1)<5 && sei(i,1)>=0)
        sei(i,1)=0;
    end
end

for i=1:1:nx
hsei(:,i) = hilbert(sei(:,1)); 
end 

% 
% x=(1:8000)*0.01;
% y=sin(x);
% hsei(:,1)=hilbert(y);

figure(1),
subplot(311),plot(1:nt,imag(hsei(:,1)),'r'),hold on 
subplot(312),plot(1:nt,real(hsei(:,1)),'b'),hold on
subplot(313),plot(1:nt,abs(atan(imag(hsei(:,1))./real(hsei(:,1)))))




