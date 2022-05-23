function ber_theory=AWGN_CE_ber(M,EbN0,LL)
%% ���룺MΪ������ EsN0���������   LLΪ2pih
%% �����ber_theoryOFDMϵͳM���Ƶ����µ����������
 ber_theory=zeros(1,length(EbN0));
 sq=M;
for m=1:length(EbN0)
    var1=(2*(sq-1))/(sq*log2(sq));
    aa=10.^(EbN0(m)/10);
    var2=sqrt((6*log2(sq)*aa)/(M^2-1));
    ber_theory(m)=var1*Q(LL*var2);
end
