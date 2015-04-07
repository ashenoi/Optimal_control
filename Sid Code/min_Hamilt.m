function u_star = min_Hamilt(P)
u_star = [];
for i = 1:1:length(P)
    if P(1,i) > 0
        u = 1;
    else
        u = 2;
    end
    u_star = [u_star;u];
end
end
