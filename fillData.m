function Yf = fillData(Y)

Yf = Y;

indmiss = find(Y<0);

aa = 0;
gap{1} = [];
prev = -1;
for jt = indmiss
    if jt-1 > prev
        aa = aa + 1;
        gap{aa} = [];
    end
    gap{aa} = [gap{aa} jt];
    prev = jt;
end


ee = 1.25;



for jg = 1:length(gap)
    if length(gap{jg}) == 1
        if gap{jg} == 1
            Yf(1) = Yf(2);
        elseif gap{jg} == length(Y)
            Yf(end) = Y(end-1);
        else
            Yf(gap{jg}) = .5*(Y(gap{jg}-1) + Y(gap{jg}+1));
        end
    else
        if ismember(length(Y),gap{jg})
            Yf(gap{jg}) = Y(min(gap{jg})-1)*ee.^(min(gap{jg})-1-gap{jg});
        elseif ismember(1,gap{jg})
            Yf(gap{jg}) = Y(max(gap{jg})+1)./ee.^(max(gap{jg})+1 - gap{jg});
        else
            Yf(gap{jg}) = Y(min(gap{jg})-1) + (gap{jg}-min(gap{jg})+1)/(length(gap{jg})+2)*(Y(max(gap{jg})+1)-Y(min(gap{jg})-1));
            Yf(gap{jg}) = Yf(gap{jg}).*(ee.^(-min(gap{jg}-min(gap{jg})+1,max(gap{jg})-gap{jg}+1)));
        end   
    end
end
        
        
        










