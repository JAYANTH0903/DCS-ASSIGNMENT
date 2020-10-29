    T = readtable('dcs.csv'); % reading the dataset

    a = T(:, 2); % Extracting the second column from the dataset which is avg temperature
    y=table2array(a)% Convarting the values which are of table type to vectors of an array

    %quabtizaiton
    f=365 % range of the dataset
    n=5 % number of bits per sample
    q=f/(2^n-1) 
    t=1:1:365 % partitioning
    x0=fix(y/q)
    y1=x0*q % quantized value
    figure(1);
    plot(t,y) 
    xlabel('Days', 'fontsize', 12, 'fontweight', 'bold')
    ylabel('Average temperature', 'fontsize', 12, 'fontweight', 'bold')
    set(gca , 'FontSize',13,'FontWeight','bold')
    title('Raw data plot','fontsize', 12, 'fontweight', 'bold')
    grid on;
    figure(2);
    plot(t,y1,"r");
    xlabel('Days', 'fontsize', 12, 'fontweight', 'bold')
    ylabel('Average temperature', 'fontsize', 12, 'fontweight', 'bold')
    set(gca , 'FontSize',13,'FontWeight','bold')
    title('Quantized data plot','fontsize', 12, 'fontweight', 'bold')
    grid on;





    x_1 = adeltamod(y1,1,1,2);%passing values to the adeltamod function

    %DELTA MODULATED SIGNAL
    figure(3);
    stairs(x_1); % to plot the modulated signal
    title('Delta modulated signal plot');
    xlabel('DAYS','fontsize',10,'fontweight','bold');
    ylabel('AVERAGE TEMPERATURE','fontsize',10,'fontweight','bold');
    grid on;

    %Demodulation /reconstruction using LOW PASS FILTER
    [val,d]=butter(1,0.2,'low');% Low pass filter function which takes order , ripple factor and type as input
    m=filter(val,d,x_1);
    figure;
    plot(m,'m');
    title('Demodulated signal with low pass filter');
    xlabel('DAYS','fontsize',10,'fontweight','bold');
    ylabel('AVERAGE TEMPERATURE','fontsize',10,'fontweight','bold');
    grid on;


    function [ADMout] = adeltamod(sig_in, Delta, td, ts)
    % Delta - minimum step size which is multiplied 2mx times when required
    % sig-in - signal input which should be a vector
    % td - orginal sampling period of the input signal - sig_in
    % ts - required sampling period of the ADM output.

        if (round(ts/td) >= 2)
            Nf = round(ts/td);    %Nearest integer
            x_sig = downsample(sig_in,Nf);
            L_sig = length(x_sig);
            L_sig_in = length(sig_in);

            ADMout = zeros(L_sig_in);    %Initialising output

            cnt1 = 0;   %Counters for no. of previous consecutively increasing 
            cnt2 = 0;   %steps
            sum = 0;
            for i=1:L_sig % to access every element of the input signal

                if (x_sig(i) == sum)
                elseif (x_sig(i) > sum)
                    if (cnt1 <=2)
                        sum = sum + Delta;  %Step up by Delta, same as in DELTA MODULATION
                    else
                        sum = sum + 2*Delta;    %Double the step size after
                                                %first two increase

                    end
                    if (sum < x_sig(i))
                        cnt1 = cnt1 + 1;
                    else
                        cnt1 = 0;
                    end
                else
                    if (cnt2 <=2)
                        sum = sum - Delta; %Step down by Delta, same as in DELTA MODULATION
                    else 
                        sum = sum - 2*Delta;


                    end
                    if (sum > x_sig(i))
                        cnt2 = cnt2 + 1;
                    else
                        cnt2 = 0;
                    end
                end
                ADMout(((i-1)*Nf + 1):(i*Nf)) = sum;
            end
        end
    end