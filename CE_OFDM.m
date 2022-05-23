%% CE-OFDMϵͳ���Ʋ�ͬM-QAM��������λ����ϵͳ��BER����
clear all
close all
clc
%% �̶������������
rng default
%% ====================================================================����
Ns=1;                                                %������
N_DFT=512;                                           %IFFT����--��������
M=64;                                                %64QAM
Norm=1/sqrt(42);                                     %64 1/sqrt(42) 16 1/sqrt(10) 4 1/sqrt(2)
bit_modu=log2(M);                                    %���ƽ�����
M_pam=sqrt(M);                                       %��Ӧ��pam���ƽ���
Tused=128e-6;                                        %��Чʱ�䳤��
N=128;                                               %���ز���                         
Nqam=N/2;                                            
Nzp=N_DFT-N-2;                                       %������
C_os=N_DFT/(N_DFT-Nzp);                              %������ϵ��
Ts=Tused/(N_DFT);                                    %�������
SNR=(0:35);                                          %�����
sample_point=round(Tused/Ts);                        %��������
MonteCarlo_num=1000;                                 %���ؿ���ʵ�����
%% ��λ���Ʋ�������
A=1;                                                 %��λ���Ƶķ�ֵ
index_2pi=1;                                         %����ϵ��2��h
phase_frist=0;                                       %��ʼ��λ
%% ===========================================================================QAM��Ϣ����                                  
mn=randi([0 1],Nqam*Ns*bit_modu,1);                  %���ɱ�����������ֱ�����ɷ��ţ���ΪҪ�Ƚϱ��ػ�������ʣ������Ϊ�����
data_convert=reshape(mn,bit_modu,Nqam*Ns);
data_convert_T=data_convert';                        %��λ����
data_decima=bi2de(data_convert_T);
a_qam=Norm*qammod(data_decima,M);
a_mn=reshape(a_qam,Nqam,Ns);%����ת��
%% Ƶ��QAM���ƺ���ź�����
sigPower_f_amn = sum(abs(a_mn(:,1)).^2)/numel(a_mn(:,1));
%% X[K]  �ź��ع�
X_K=zeros(N_DFT,Ns);
for m=1:Ns
    one=a_mn(:,m);
    two=zeros(length(one),1);
    two_index=fliplr(1:length(one));%%������ͷβת��
    for i=1:length(one)
        two(i)=one(two_index(i));
    end
    two=conj(two);
    X_K(:,m)=[0;one;zeros(Nzp,1);0;two];
end
%% Ƶ��X[K]�ź�����
sigPower_f_xk = sum(abs(X_K(:,1)).^2)/numel(X_K(:,1));
%% ======================================================== IDFT>>����m(t)==================����m(t)
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
%% m(t)  20220307 д������
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
%% =========================================================== IDFT>>>����x[n]================���Ĺ�ʽ��1��
%% xt
xt=zeros(N_DFT,Ns);
IDFT_nor_xt=(1/(N_DFT))*2;
% �����ź�ģ��
for m=0:Ns-1
    for t_point=0:N_DFT-1
        temp=0;
        for k=1:Nqam
            temp=temp+(1/(N_DFT))*2*(real(X_K(k+1,m+1))*cos(2*pi*k*t_point/(N_DFT))-imag(X_K(k+1,m+1))*sin(2*pi*k*t_point/(N_DFT)));
        end
        xt(t_point+1,m+1)=temp;
    end
end
%% xt�ź�����
sigPower_xt = sum(abs(xt(:,1)).^2)/numel(xt(:,1));
%% xt��һ����������
xt=xt/sqrt(sigPower_xt);
sigPower_xt_norm = sum(abs(xt(:,1)).^2)/numel(xt(:,1));
%% ======================================================Phase modulation
st=zeros(sample_point,Ns);%�ռ�Ԥ����
for m=0:Ns-1
    for i=0:sample_point-1
        st(i+1,m+1)=A*exp(1j*(index_2pi*(real(xt(i+1,m+1))))+phase_frist);%%������ǰ�����ɵ�m_t��x[n]ѡ������һ��
    end
end
%% ================================================================�ŵ��ͽ���=======================
BER_get=zeros(1,length(SNR)); 
MontCarspace=zeros(1,MonteCarlo_num);
for k=1:length(SNR)
    for MC=1:MonteCarlo_num
        r_t=awgn(st,SNR(k));
        %% ���ն��ź����� ����������
%         SNR=20;
%         % �źŹ��� ȡ��һ�����ŵĲ�������
%         sigPower = sum(abs(st(:,1)).^2)/numel(st(:,1));
%         % ��������
%         noise_power=sigPower/(10^(SNR/10));
        %% ������ �����
        ss=zeros(sample_point,Ns);
        for i=1:Ns
            cccc=phase(r_t(:,i),pi);
            ss(:,i)=unwrap(cccc)/(index_2pi);
        end
        %% ���ն˽���λ�������ź�ȥ��һ��
        ss=ss*sqrt(sigPower_xt);
        %% ���ն�OFDM���
        xt_f=fft(ss);%%  �Ż���������д����
        %% CE�źŽ⹹
        X_R_dec=zeros(Nqam,Ns);
        for i=1:Ns
            X_R_dec(:,i)=xt_f(2:(Nqam+1),i);
        end
        %% QAM���
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

legend('CE-OFDM','CE-OFDM����ֵ')

xxxx=1;



