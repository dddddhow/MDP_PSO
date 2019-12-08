clc,clear;
%����������ݣ�sei����ȡ����������
nt=8000;
nx=200;
dt=0.0004;
sei=zeros(nt,nx);
fp=fopen('record','rb');
sei=fread(fp,[nt nx],'float32');
fclose(fp);
% imagesc(sei)

% %һ������˲����(�������amp����λ����ang��Ƶ������fre)
% hsei=zeros(nt,nx);
% hsei=hilbert(sei);
% amp=sqrt(real(hsei).*real(hsei)+imag(hsei).*imag(hsei));
% figure(1)
% imagesc(amp),xlabel('����/��'),ylabel('ʱ��'),title('�������');
% ang=atan(imag(hsei)./real(hsei));
% figure(2)
% imagesc(ang),xlabel('����/��'),ylabel('ʱ��'),title('��λ����');
% fre=diff(ang,1)/(2*pi);
% figure(3)
% imagesc(fre),xlabel('����/��'),ylabel('ʱ��'),title('Ƶ������');


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




