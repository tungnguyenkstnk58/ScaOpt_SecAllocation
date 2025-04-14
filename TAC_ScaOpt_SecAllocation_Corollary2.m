close all
clear
clc
%%
expN = 1; % the number of experiments
exptimeCYC = zeros(expN, 1); % computation time by Theorem 1
exptimeKYP = zeros(expN, 1); % computation time by Corollary 2
expcostCYC = zeros(expN, 1); % worst-case disruption by Theorem 1
expcostKYP = zeros(expN, 1); % worst-case disruption by Corollary 2
%% start experiment
for expk = 1 : expN
    rng(10 + expk); % choose random seed
    %% System parameters
    N = 50; % network size
    [L, A] = ConnectedGraph(N,0); % Laplacian matrix L and adjacency matrix A
    Theta = 1.5 * eye(N); % self-loop gains
    L = L + Theta; % weighted Laplacian with self-loops
    w = 2 * ones(N,1); % performance weighting factor
    beta = 3; % sensor budget
    delta = 0.7 * ones(N,1); % alarm threshold
    kappa = 0.1 * ones(N,1); % sensor cost
    na = 3; % attack types
    alpha = 1 : na; % attack budget vector
    phi = [0.5, 0.35, 0.15]; % attack type probability 
    E = 1e1; % maximum attack energy
    nAs = 0; % total numbers attack scenarios
    for k = 1:na
        nAs = nAs + nchoosek(N, alpha(k)); % attack scenario for each attack type
    end

    %% Security Allocation Solution - Theorem 1
    Mt = 1e5; % big M
    % all the optimal variables in Theorem 1
    zM = binvar(N, 1); % optimal variable for monitor set M
    omegaM = cell(nAs, 1); 
    P = cell(nAs, 1); 
    Qk = cell(na, 1);
    psi = cell(nAs, 1);
    QAn = cell(nAs, 1); 
    BAn = cell(nAs, 1);
    As_pre = 0;

    % constraints
    F = [];
    F = [F, ones(1,N) * zM <= beta]; % sensor budget
    % objective 
    h = 0;

    for k = 1 : na
        Qk{k} = sdpvar(1, 1);
        F = [F, Qk{k} >= 0];
        As_sum = nchoosek(1 : N, alpha(k)); % all the scenarios for attack type k
        for As = 1 : size(As_sum,1) % for each attack set of attack type k
            An = As_sum(As,:); % take an attack set
            BAn{As_pre+As} = [];
            for ka = 1 : size(An, 2)
                ea = zeros(N, 1);
                ea(An(ka)) = 1;
                BAn{As_pre+As} = [BAn{As_pre+As}, ea]; % create attack entry matrix based on the taken attack set
            end
            % optimal variable for each attack set
            omegaM{As_pre+As} = sdpvar(N, 1);
            psi{As_pre+As} = sdpvar(size(An,2), 1);
            P{As_pre+As} = sdpvar(N, N, 'symmetric');
    
            % constraints in Theorem 1
            F = [F, omegaM{As_pre+As} >= 0];
            F = [F, psi{As_pre+As} >= 0];
            F = [F, omegaM{As_pre+As} <= Mt * zM];
            QAn{As_pre+As} = delta' * omegaM{As_pre+As} + E * ones(1,size(An,2)) * psi{As_pre+As};
            F = [F, QAn{As_pre+As} <= Qk{k}];
            lmi = [-L'*P{As_pre+As}-P{As_pre+As}*L P{As_pre+As}*BAn{As_pre+As};
                   BAn{As_pre+As}'*P{As_pre+As} -diag(psi{As_pre+As})] + diag([w.^2; zeros(size(An,2),1)]) - diag([omegaM{As_pre+As}; zeros(size(An,2),1)]);
            F = [F, lmi <= 0];
        end
        h = h + phi(k) * (kappa' * zM + Qk{k});
        As_pre = As_pre + As;
    end

% optimization
options = sdpsettings('verbose', 2, 'debug', 0, 'solver', 'bnb', 'bnb.solver', 'mosek');
CYC = optimize(F, h, options);
exptimeCYC(expk) = CYC.solvertime;
expcostCYC(expk) = value(h);

%% Security Allocation Solution by Corollary 2
% constraints
F = [];
F = [F, ones(1,N) * zM <= beta;]; % sensor budget
% objective 
h = 0;

for k = 1 : na
    Qk{k} = sdpvar(1, 1);
    F = [F, Qk{k} >= 0];
    As_sum = nchoosek(1 : N, alpha(k)); % all the scenarios for attack type k
    for As = 1 : size(As_sum,1) % for each attack set of attack type k
        An = As_sum(As,:); % take an attack set
        BAn{As_pre+As} = [];
        for ka = 1 : size(An,2)
            ea = zeros(N, 1);
            ea(An(ka)) = 1;
            BAn{As_pre+As} = [BAn{As_pre+As}, ea]; % create attack entry matrix based on the taken attack set
        end
        % optimal variable for each attack set
        omegaM{As_pre+As} = sdpvar(N, 1);
        psi{As_pre+As} = sdpvar(size(An,2), 1);
        P{As_pre+As} = diag(sdpvar(N, 1)); % diagonal P

        % constraints
        F = [F, omegaM{As_pre+As} >= 0];
        F = [F, psi{As_pre+As} >= 0];
        F = [F, omegaM{As_pre+As} <= Mt * zM];
        F = [F, P{As_pre+As} >= 0];
        QAn{As_pre+As} = delta' * omegaM{As_pre+As} + E * ones(1,size(An,2)) * psi{As_pre+As};
        F = [F, QAn{As_pre+As} <= Qk{k}];
        lmi = [-L'*P{As_pre+As}-P{As_pre+As}*L P{As_pre+As}*BAn{As_pre+As};
               BAn{As_pre+As}'*P{As_pre+As} -diag(psi{As_pre+As})] + diag([w.^2; zeros(size(An,2),1)]) - diag([omegaM{As_pre+As}; zeros(size(An,2),1)]);
        F = [F, lmi <= 0];
    end
    h = h + phi(k) * (kappa' * zM + Qk{k});
    As_pre = As_pre + As;
end

% optimization
KYP = optimize(F,h,options);
exptimeKYP(expk) = KYP.solvertime;
expcostKYP(expk) = value(h);

end

figure(1)
boxplot(log(expcostKYP./expcostCYC),'Labels',{'N = 10'},'whisker',inf);
ylabel('$$\log(R^\star_{T2}/R^\star_{T1})$$','Interpreter','latex','FontSize',14);
title('Optimal defense cost comparison','FontSize',14);
set(gca,'FontSize',14,'fontWeight','bold');
%print(gcf,'figs/10NodeKYPcomparecost','-dpng','-r300');

figure(2)
boxplot(exptimeKYP./exptimeCYC*100,'Labels',{'N = 10'},'whisker',inf);
ylabel('$$T_2/T_1(\%)$$','Interpreter','latex','FontSize',14);
title('Computational cost comparison','FontSize',14);
set(gca,'FontSize',14,'fontWeight','bold');
%print(gcf,'figs/10NodeKYPcomparetime','-dpng','-r300');
% plot(exptimeCYC,'--','LineWidth',2);
% hold on
% plot(exptimeKYP,'-.','LineWidth',2);
% set(gca,'Fontsize',16);
% ylabel('Computational time (s)','FontSize',16);
% xlabel('Experiment','FontSize',16);
% xticks(0:5:expN);
% xlim([1 expN]);
% ylim([0 max(exptimeCYC)*1.3]);
% grid;
% legend('Cyclo dissipative','KYP');
%print(gcf,'figs/10NodeKYPcompare','-dpng','-r300');


%% Test
% omegaM_max = 0;
% for k = 1:size(Akv,1)
%     if max(value(omegaM{k})) > omegaM_max
%         omegaM_max = max(value(omegaM{k}));
%         k
%     end
% end




