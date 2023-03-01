clc
clear 
close all

%% System model of received signal
K=1000;k=2; m=4;
u = randi([0 1],1,K); % generate a sequence of K random bits, (0/1)

v=encoder(u,2); %encoder.m that takes a sequence u as input and maps it to a sequence v
v_prime=randintrlv(v,rand(1,1));%The next step is to randomly scramble (a.k.a. interleave or permute) the encoded sequence v
N=2*(K+3);
group=N/k; % In total we must have N/log2M groups.
%% Gray Mapping

a_gray= graymapping(v_prime,k,N);

%Adding Noise
EbN0dB=-5;
    error_hard_graymapping=zeros(16,1);
    error_hard_nongray=zeros(16,1);
    error_soft_gray=zeros(16,1);
    error_soft_nongray=zeros(16,1);
while(EbN0dB<=10)
    EbN0=10^(EbN0dB/10); 
    w=sqrt(1/(EbN0))*randn(1,length(a_gray));
    r_gray=a_gray+w';   
    
    %Lk Sequence Calculatin
    Lk=mydemodulator(r_gray,k,N);
    Lk=randdeintrlv(Lk,rand(1,1)); % inverse scrambler
    Lk=Lk';
    udot_gray=out_decoder(2,Lk); %takes a sequence of values Lk and finds the most likely sequence udot
    
    %% nongraymapping
    a_nongraymapping= nongraymapping(v_prime,k,N);
    %Adding Noise
    r_nongray=a_nongraymapping+w';   
    Lk_nongray=mydemodnongray(r_nongray,k,N);
    Lk_nongray=randdeintrlv(Lk_nongray,rand(1,1));
    Lk_nongray=Lk_nongray'; % inverse scrambler
    udot_nongray=out_decoder(2,Lk_nongray); %takes a sequence of values Lk and finds the most likely sequence udot
    %soft part
    udot_soft_gray= soft(r_gray,N);
    udot_soft_nongray= soft_nongray(r_nongray,N);
    %error calculation
    error_hard_graymapping(EbN0dB+6) = errorcalculation(udot_gray,u,K);
    error_hard_nongray(EbN0dB+6)=errorcalculation(udot_nongray,u,K);
    error_soft_gray(EbN0dB+6)=errorcalculation(udot_soft_gray,u,K);
    error_soft_nongray(EbN0dB+6)=errorcalculation(udot_soft_nongray,u,K);
    EbN0dB=EbN0dB+1;
  
end
EbN0dB_axis=-5:1:10;
 figure
 semilogy(EbN0dB_axis,error_hard_graymapping);
 hold on
 semilogy(EbN0dB_axis,error_hard_nongray);
  title('BER curves versus Eb/N0');
 xlabel('Eb/NO(dB)');
 ylabel('BER');
legend('hardgraymapping)','hardnongray')
figure
semilogy(EbN0dB_axis,error_soft_gray);
hold on
 semilogy(EbN0dB_axis,error_soft_nongray);
 title('BER curves versus Eb/N0');
 xlabel('Eb/NO(dB)');
 ylabel('BER');
legend('softgray','softnongray');
%%figure all

figure
semilogy(EbN0dB_axis,error_hard_graymapping);
 hold on
 semilogy(EbN0dB_axis,error_hard_nongray);
 hold on
semilogy(EbN0dB_axis,error_soft_gray);
hold on
 semilogy(EbN0dB_axis,error_soft_nongray);
 title('BER curves versus Eb/N0');
 xlabel('Eb/NO(dB)');
 ylabel('BER');
legend('hardgray','hardnongray','softgray','softnongray');

function a_gray= graymapping(v_prime,k,N)
i=1; j=1;
group=N/k; % In total we must have N/log2M groups.
a_gray=zeros(group,1);
while(i<=length(v_prime))
    if(v_prime(i)==0 && v_prime(i+1)==0)
        a_gray(j)=-3;
    end
    if(v_prime(i)==0 && v_prime(i+1)==1)
        a_gray(j)=-1;
    end
    if(v_prime(i)==1 && v_prime(i+1)==1)
        a_gray(j)=1;
    end
    if(v_prime(i)==1 && v_prime(i+1)==0)
        a_gray(j)=3;
    end
    i=i+2;
    j=j+1;
end
end
function a_nongray = nongraymapping(v_prime,k,N)
i=1; j=1;
group=N/k; % In total we must have N/log2M groups.
a_nongray=zeros(group,1);
while(i<=length(v_prime))
    if(v_prime(i)==0 && v_prime(i+1)==0)
        a_nongray(j)=-3;
    end
    if(v_prime(i)==0 && v_prime(i+1)==1)
        a_nongray(j)=-1;
    end
    if(v_prime(i)==1 && v_prime(i+1)==0)
        a_nongray(j)=1;
    end
    if(v_prime(i)==1 && v_prime(i+1)==1)
        a_nongray(j)=3;
    end
    i=i+2;
    j=j+1;
end
end
function Lk = mydemodulator(receivedsignal,k,N) 


E0=9; E1=1; E2=1; E3=9;
z0=-3; z1=-1; z2=1; z3=3;
group=N/k;
i=1;
temp=zeros(4,1);
temp_1=zeros(N/2,1);
while(i<=group) %%dedection of location of the decisions
    temp(1)=receivedsignal(i)*z0-(E0/2);
    temp(2)=receivedsignal(i)*z1-(E1/2);
    temp(3)=receivedsignal(i)*z2-(E2/2);
    temp(4)=receivedsignal(i)*z3-(E3/2);
    [~,I] = max(temp);
    temp_1(i)=I-1;
    i=i+1;

end

i=1;
temp_2=zeros(length(temp_1),1);
while(i<=group) %% dedection of decision variables
    if(temp_1(i)==0)
        temp_2(i)=-3;
    end
    if(temp_1(i)==1)
        temp_2(i)=-1;
    end
    if(temp_1(i)==2)
        temp_2(i)=1;
    end
    if(temp_1(i)==3)
        temp_2(i)=3;
    end
    i=i+1;
end

i=1; j=1;
vdotprime=zeros(N,1);
while(i<=group) %%demapping of decision variables
    if(temp_1(i)==0)
        vdotprime(j)=0;
        vdotprime(j+1)=0;
    end
    if(temp_1(i)==1)
        vdotprime(j)=0;
        vdotprime(j+1)=1;
    end
    if(temp_1(i)==2)
        vdotprime(j)=1;
        vdotprime(j+1)=1;
    end
    if(temp_1(i)==3)
        vdotprime(j)=1;
        vdotprime(j+1)=0;
    end
    i=i+1;
    j=j+2;
end

i=1;
Lk=zeros(N,1);
while(i<=N) %%calculation of L
    if(vdotprime(i)==1)
        Lk(i)=10;
    end
    if(vdotprime(i)==0) %temp2=demoduladion result before demapping, r=received signal
        Lk(i)=-10;
    end
    i=i+1;
end
end

function lknongray = mydemodnongray(receivedsignal,k,N)


E0=9; E1=1; E2=1; E3=9;
z0=-3; z1=-1; z2=1; z3=3;
group=N/k;
i=1;
temp=zeros(4,1);
temp_1=zeros(N/2,1);
while(i<=group) %%dedection of location of the decisions
    temp(1)=receivedsignal(i)*z0-(E0/2);
    temp(2)=receivedsignal(i)*z1-(E1/2);
    temp(3)=receivedsignal(i)*z2-(E2/2);
    temp(4)=receivedsignal(i)*z3-(E3/2);
    [~,I] = max(temp);
    temp_1(i)=I-1;
    i=i+1;

end

i=1;
temp_2=zeros(length(temp_1),1);
while(i<=group) %% dedection of decision variables
    if(temp_1(i)==0)
        temp_2(i)=-3;
    end
    if(temp_1(i)==1)
        temp_2(i)=-1;
    end
    if(temp_1(i)==2)
        temp_2(i)=1;
    end
    if(temp_1(i)==3)
        temp_2(i)=3;
    end
    i=i+1;
end

i=1; j=1;
vdotprime=zeros(N,1);
while(i<=group) %%demapping of decision variables
    if(temp_1(i)==0)
        vdotprime(j)=0;
        vdotprime(j+1)=0;
    end
    if(temp_1(i)==1)
        vdotprime(j)=0;
        vdotprime(j+1)=1;
    end
    if(temp_1(i)==2)
        vdotprime(j)=1;
        vdotprime(j+1)=0;
    end
    if(temp_1(i)==3)
        vdotprime(j)=1;
        vdotprime(j+1)=1;
    end
    i=i+1;
    j=j+2;
end

i=1;
lknongray=zeros(N,1);
while(i<=N) %%calculation of L
    if(vdotprime(i)==1)
        lknongray(i)=10;
    end
    if(vdotprime(i)==0) %temp2=demoduladion result before demapping, r=received signal
        lknongray(i)=-10;
    end
    i=i+1;
end
end

function udot_soft = soft(receivedsignal,N) 


i=1;j=1;
Lprime=zeros(N,1);
while(i<=length(Lprime))

    temp1= exp(-1*((receivedsignal(j)-(-3))^2)) + exp(-1*((receivedsignal(j)-(-1))^2)); %zero position at 0 and 1
    temp2= exp(-1*((receivedsignal(j)-(1))^2)) + exp(-1*((receivedsignal(j)-(3))^2));%one position at 2 and 3

    Lprime(i)= log(temp2/temp1);
    i=i+1;

    temp1= exp(-1*((receivedsignal(j)-(-3))^2)) + exp(-1*((receivedsignal(j)-(3))^2));
    temp2= exp(-1*((receivedsignal(j)-(-1))^2)) + exp(-1*((receivedsignal(j)-(1))^2));
    Lprime(i)= log(temp2/temp1);
    i=i+1;
    j=j+1;
end
Lprime_new=randdeintrlv(Lprime,rand(1,1));
Lprime_new=Lprime_new'; % inverse scrambler
udot_soft=out_decoder(2,Lprime_new);

end

function udot_soft_nongray = soft_nongray(receivedsignal,N) 
i=1;j=1;
Lprime=zeros(N,1);
while(i<=length(Lprime))

    temp1= exp(-1*((receivedsignal(j)-(-3))^2)) + exp(-1*((receivedsignal(j)-(-1))^2)); %zero position at 0 and 1
    temp2= exp(-1*((receivedsignal(j)-(1))^2)) + exp(-1*((receivedsignal(j)-(3))^2));%one position at 2 and 3

    Lprime(i)= log(temp2/temp1);
    i=i+1;

    temp1= exp(-1*((receivedsignal(j)-(-3))^2)) + exp(-1*((receivedsignal(j)-(1))^2));
    temp2= exp(-1*((receivedsignal(j)-(-1))^2)) + exp(-1*((receivedsignal(j)-(3))^2));
    Lprime(i)= log(temp2/temp1);
    i=i+1;
    j=j+1;
end
Lprime_new=randdeintrlv(Lprime,rand(1,1));
Lprime_new=Lprime_new'; % inverse scrambler
udot_soft_nongray=out_decoder(2,Lprime_new);

end
function error = errorcalculation(udot,u,K)


error=0;i=1;
while(i<=K)
    if(udot(i)~=u(i))
        error=error+1;
    end
    i=i+1;
end
error=(error/K)*100;
end