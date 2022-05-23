%% CE-OFDM系统绘制不同M-QAM调制下相位调制系统的BER曲线
clear all
close all
clc
%% 固定产生的随机数
rng default
%% ====================================================================参数
Ns=1;                                                %符号数
N_DFT=512;                                           %IFFT点数--采样点数
M=64;                                                %64QAM
Norm=1/sqrt(42);                                     %64 1/sqrt(42) 16 1/sqrt(10) 4 1/sqrt(2)
bit_modu=log2(M);                                    %调制进制数
M_pam=sqrt(M);                                       %对应的pam调制阶数
Tused=128e-6;                                        %有效时间长度
N=128;                                               %子载波数                         
Nqam=N/2;                                            
Nzp=N_DFT-N-2;                                       %补零数
C_os=N_DFT/(N_DFT-Nzp);                              %过采样系数
Ts=Tused/(N_DFT);                                    %采样间隔
SNR=(0:35);                                          %信噪比
sample_point=round(Tused/Ts);                        %采样点数
MonteCarlo_num=1000;                                 %蒙特卡洛实验次数
%% 相位调制参数设置
A=1;                                                 %相位调制的幅值
index_2pi=1;                                         %调制系数2Πh
phase_frist=0;                                       %初始相位
%% ===========================================================================QAM信息产生                                  
mn=randi([0 1],Nqam*Ns*bit_modu,1);                  %生成比特流而不是直接生成符号，因为要比较比特获得误码率，解调端为逆过程
data_convert=reshape(mn,bit_modu,Nqam*Ns);
data_convert_T=data_convert';                        %高位在右
data_decima=bi2de(data_convert_T);
a_qam=Norm*qammod(data_decima,M);
a_mn=reshape(a_qam,Nqam,Ns);%串并转换
%% 频域QAM调制后的信号能量
sigPower_f_amn = sum(abs(a_mn(:,1)).^2)/numel(a_mn(:,1));
%% X[K]  信号重构
X_K=zeros(N_DFT,Ns);
for m=1:Ns
    one=a_mn(:,m);
    two=zeros(length(one),1);
    two_index=fliplr(1:length(one));%%将序列头尾转换
    for i=1:length(one)
        two(i)=one(two_index(i));
    end
    two=conj(two);
    X_K(:,m)=[0;one;zeros(Nzp,1);0;two];
end
%% 频域X[K]信号能量
sigPower_f_xk = sum(abs(X_K(:,1)).^2)/numel(X_K(:,1));
%% ======================================================== IDFT>>生成m(t)==================论文m(t)
%% I_K
I_K=zeros(N,Ns);
for m=1:Ns
    for k=1:N
        if k<=Nqam
            I_K(k,m)=real(a_mn(k,m));
        else
            I_K(k,m)=-imag(a_mn(k-Nqam,m));
        end
    end
end
%% m(t)  20220307 写到这里
IDFT_nor_mt=(Ts/Tused)*2;
m_t=zeros(N_DFT,Ns);
for m=0:Ns-1
    for t_index=0:N_DFT-1
        temp=0;
        for k=1:N
            if k<=Nqam
                temp=temp+IDFT_nor_mt*I_K(k,m+1)*cos(2*pi*k*t_index*Ts/(Tused));%
            else
                temp=temp+IDFT_nor_mt*I_K(k,m+1)*sin(2*pi*(k-Nqam)*t_index*Ts/Tused);
            end
        end
        m_t(t_index+1,m+1)=temp;
    end
end
%% =========================================================== IDFT>>>生成x[n]================论文公式（1）
%% xt
xt=zeros(N_DFT,Ns);
IDFT_nor_xt=(1/(N_DFT))*2;
% 建立信号模型
for m=0:Ns-1
    for t_point=0:N_DFT-1
        temp=0;
        for k=1:Nqam
            temp=temp+(1/(N_DFT))*2*(real(X_K(k+1,m+1))*cos(2*pi*k*t_point/(N_DFT))-imag(X_K(k+1,m+1))*sin(2*pi*k*t_point/(N_DFT)));
        end
        xt(t_point+1,m+1)=temp;
    end
end
%% xt信号能量
sigPower_xt = sum(abs(xt(:,1)).^2)/numel(xt(:,1));
%% xt归一化后信能量
xt=xt/sqrt(sigPower_xt);
sigPower_xt_norm = sum(abs(xt(:,1)).^2)/numel(xt(:,1));
%% ======================================================Phase modulation
st=zeros(sample_point,Ns);%空间预分配
for m=0:Ns-1
    for i=0:sample_point-1
        st(i+1,m+1)=A*exp(1j*(index_2pi*(real(xt(i+1,m+1))))+phase_frist);%%在这里前面生成的m_t和x[n]选择其中一种
    end
end
%% ================================================================信道和接收=======================
BER_get=zeros(1,length(SNR)); 
MontCarspace=zeros(1,MonteCarlo_num);
for k=1:length(SNR)
    for MC=1:MonteCarlo_num
        r_t=awgn(st,SNR(k));
        %% 接收端信号能量 与噪声能量
%         SNR=20;
%         % 信号功率 取第一个符号的采样点数
%         sigPower = sum(abs(st(:,1)).^2)/numel(st(:,1));
%         % 噪声功率
%         noise_power=sigPower/(10^(SNR/10));
        %% 反正切 解卷绕
        ss=zeros(sample_point,Ns);
        for i=1:Ns
            cccc=phase(r_t(:,i),pi);
            ss(:,i)=unwrap(cccc)/(index_2pi);
        end
        %% 接收端将相位解调后的信号去归一化
        ss=ss*sqrt(sigPower_xt);
        %% 接收端OFDM解调
        xt_f=fft(ss);%%  优化，后续自写函数
        %% CE信号解构
        X_R_dec=zeros(Nqam,Ns);
        for i=1:Ns
            X_R_dec(:,i)=xt_f(2:(Nqam+1),i);
        end
        %% QAM解调
        pro_demodu_qam=reshape(X_R_dec,1,Nqam*Ns);
        demodu_qam=qamdemod(pro_demodu_qam/Norm,M);
        demodu_qam_convert=de2bi(demodu_qam);
        demo_data_convert_T=demodu_qam_convert';
        mu_demodu=reshape(demo_data_convert_T,Nqam*Ns*bit_modu,1);
        [val,ber]=biterr(mu_demodu,mn);
        MontCarspace(1,MC)=ber;
    end
    ber_MC=mean(MontCarspace);
    BER_get(1,k)=ber_MC;
end
%% CE-OFDM result of simulation
terory=AWGN_CE_ber(sqrt(M),SNR,index_2pi);
close all
figure(201) 
semilogy(SNR,BER_get,'-*')
ylim([1e-4 1])
xlim([0 35])
title('Bit error ratio(M-QAM)')
xlabel('SNR')
ylabel('BER')
grid on;
hold on
semilogy(SNR,terory,'-o');

legend('CE-OFDM','CE-OFDM理论值')

xxxx=1;



