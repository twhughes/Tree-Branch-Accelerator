tyler_study;


NT = 100;

[~,L_list] = s.get_Length_list(1);
L_list=  L_list(2:end);
NL = length(L_list);

tau_list = 10.^linspace(-14,-9,NT);


G_list = zeros(NL,1);

for li = (1:NL)
    s.L = L_list(li);
    s.solve_geometry;
    Q_tau;  
    G_list(li) = max(max(abs(Gs)));    
end
