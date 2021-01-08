% Iterative water-filling (TPC+PAC)
% capacity under the joint TPC+PAC constraints for a massive MIMO channel under
% favorable propagation
% inputs: 
% W : channel gains
% PT : total transmit power constraint
% P1 : per-antenna power constraints


function iterative_water_filling(W,PT,P1)
if PT >= sum(P1)
    final_power=P1;
else
    ff=1;
    ttt=2;
    py=1;
    while(ttt>1)
        k=1;
        m=length(W);
        bisectionerror=1;
        % To find mu, we can use an analytical solution or the bisection algorithm.
        % Here, the bisection algorithm is used.
        xl=0;
        xu=max(W);
        f=0;
        clear r
        while(bisectionerror>1e-8)
            midpoint=(1/2)*(xl+xu);
            for k=1:m
                a=max(0,(midpoint)^(-1)-(W(k)^(-1)));
                r(k)=a;
            end
            f_value=sum(r)-PT;
            if f_value<0
                xu=midpoint;
            elseif f_value>0
                xl=midpoint;
            end
            bisectionerror=abs(f_value);
        end

        % set of streams that exceed the PACs
        b=0;
        pac=0;
        sdsd=1;
        clear delete_W
        for i=1:m
            if r(i)>P1(i)
                r(i)=P1(i);
                pac=1+pac;
                final_W(py)=W(i);
                final_power(py)=r(i);
                delete_W(sdsd)=i;
                sdsd=1+sdsd;
                py=py+1;
            end
        end
        
        if pac>0
            pac_completed=sum(P1(delete_W));
            W(delete_W)=[];
            P1(delete_W)=[];
        end

        % pac is the number of elements in the set of streams that exceed the PACs
        if pac==0
            ttt=0;
            u_EWF=midpoint;
            for i=1:m
                final_W(py)=W(i);
                final_power(py)=r(i);
                py=py+1;
            end
        end
        
        if (ff==1 && pac==0)
            ttt=0;
        else
            PT=PT-pac_completed;
        end
        ff=ff+1;
    end
    end
final_W;
final_power;
save('PAC_TPC_EWF.mat')
end
