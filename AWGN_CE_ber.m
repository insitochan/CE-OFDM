function ber_theory=AWGN_CE_ber(M,EbN0,LL)
%% 输入：M为进制数 EsN0符号信噪比   LL为2pih
%% 输出：ber_theoryOFDM系统M进制调制下的理论信噪比
 ber_theory=zeros(1,length(EbN0));
 sq=M;
for m=1:length(EbN0)
    var1=(2*(sq-1))/(sq*log2(sq));
    aa=10.^(EbN0(m)/10);
    var2=sqrt((6*log2(sq)*aa)/(M^2-1));
    ber_theory(m)=var1*Q(LL*var2);
end
