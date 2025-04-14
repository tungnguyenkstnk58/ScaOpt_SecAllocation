close all
clear
clc
%%
figure(1) % create a figure
expN = 1; % the number of experiments
exptimeGS = zeros(expN, 1); % computation time for using Theorem 1 
exptimeBF = zeros(expN, 1); % computation time for using Brute-force
expRd = cell(expN, 1); % defense cost by two methods
%% start experiment
for expk = 1 : expN
    rng(expk + 10); % choose random seed
    expRd{expk} = []; % defense cost initializtion
    %% System parameters
    N = 10; % network size
    [L, A] = ConnectedGraph(N, 1); % Laplacian matrix L and adjacency matrix A
    Theta = 0.7 * eye(N); % self-loop gains
    L = L + Theta; % weighted Laplacian with self-loops
    w = ones(N, 1); % performance weighting factor
    beta = 3; % sensor budget
    delta = 0.5 * ones(N, 1); % alarm threshold
    kappa = 0.3 * ones(N, 1); % sensor cost
    na = 3; % attack types
    alpha = 1 : na; % attack budget vector
    phi = [0.5, 0.35, 0.15]; % attack type probability 
    E = 1e1; % maximum attack energy
    nAs = 0; % total numbers attack scenarios
    for ka = 1 : na
        nAs = nAs + nchoosek(N, alpha(ka)); % attack scenario for each attack type
    end

    %% Computing Vinf_hat by using Lemma 4
    % we first create an optimizer to compute each V_inf(A)
    psiv = sdpvar(alpha(end), 1); 
    Pav = sdpvar(N, N, 'symmetric');
    Bav = sdpvar(N, alpha(end));
    F = [];
    F = [F, psiv >= 0];
    lmi = [-L'*Pav-Pav*L+diag(w.^2) Pav*Bav;
            Bav'*Pav -diag(psiv)];
    F = [F, lmi <= 0];
    h = E * ones(1, alpha(end)) * psiv;
    options = sdpsettings('verbose', 0, 'debug', 0, 'solver', 'mosek');
    Vinf = optimizer(F, h, options, Bav, h);
    Akv = nchoosek(1 : N, alpha(end));
    Vinf_max = 0;
    for k = 1 : size(Akv,1)
        Ba = zeros(N, alpha(end));
        for m = 1 : alpha(end)
            Ba(Akv(k,m), m) = 1;
        end
        Vinf_value = Vinf(Ba);
        if Vinf_value > Vinf_max
            Vinf_max = Vinf_value;
        end
    end
    Vinf_hat = Vinf_max/min(delta); % Vinf_hat after using Lemma 4

    %% Security Allocation Solution - Theorem 1
    % all the optimal variables in Theorem 1
    zM = binvar(N, 1); % optimal variable for monitor set M
    omegaM = cell(nAs,1); 
    P = cell(nAs,1); 
    Qk = cell(na,1);
    psi = cell(nAs,1);
    QAn = cell(nAs,1); 
    BAn = cell(nAs,1);
    As_pre = 0;

    % constraints
    F = [];
    F = [F, ones(1,N) * zM <= beta;]; % sensor budget
    % objective 
    h = 0;
    for k = 1 : na
        Qk{k} = sdpvar(1,1);
        F = [F, Qk{k} >= 0];
        As_sum = nchoosek(1 : N, alpha(k)); % all the scenarios for attack type k
        for As = 1 : size(As_sum,1) % for each attack set of attack type k
            An = As_sum(As,:); % take an attack set
            BAn{As_pre+As} = [];
            for ka = 1 : size(An, 2)
                ea = zeros(N,1);
                ea(An(ka)) = 1;
                BAn{As_pre+As} = [BAn{As_pre+As}, ea]; % create attack entry matrix based on the taken attack set
            end
            % optimal variable for each attack set
            omegaM{As_pre+As} = sdpvar(N, 1);
            psi{As_pre+As} = sdpvar(size(An, 2), 1);
            P{As_pre+As} = sdpvar(N, N, 'symmetric');

            % constraints in Theorem 1
            F = [F, omegaM{As_pre+As} >= 0];
            F = [F, psi{As_pre+As} >= 0];
            F = [F, omegaM{As_pre+As} <= Vinf_hat * zM];
            QAn{As_pre+As} = delta' * omegaM{As_pre+As} + E * ones(1,size(An,2)) * psi{As_pre+As};
            F = [F, QAn{As_pre+As} <= Qk{k}];
            lmi = [-L'*P{As_pre+As}-P{As_pre+As}*L, P{As_pre+As}*BAn{As_pre+As};
                    BAn{As_pre+As}'*P{As_pre+As}, -diag(psi{As_pre+As})] + diag([w.^2; zeros(size(An,2),1)]) - diag([omegaM{As_pre+As}; zeros(size(An,2),1)]);
            F = [F, lmi <= 0];
        end
        h = h + phi(k) * (kappa' * zM + Qk{k}); % objective function in Theorem 1
        As_pre = As_pre + As; % counter for the number of attack sets
    end

    % solving the MISDP in Theorem 1
    options = sdpsettings('verbose', 2, 'debug', 0, 'solver', 'bnb');
    options.bnb.solver = 'mosek';
    options.mosek.MSK_IPAR_NUM_THREADS = 8;
    sol = optimize(F, h, options);
    exptimeGS(expk) = sol.solvertime;
    value(h)
    % record experimental results
    expRd{expk} = [expRd{expk}; value(h)];
    plot(expk, expRd{expk}(1), '*', 'MarkerEdgeColor', 'red', 'MarkerSize', 16);
    hold on

    %% Security Allocation Solution by Brute-Force
    % optimal variables 
    omegaM = sdpvar(N, 1);
    P = sdpvar(N, N, 'symmetric'); 
    zMv = sdpvar(N, 1);   
    options = sdpsettings('verbose', 0, 'debug', 0, 'solver', 'mosek');
    R_matrix = cell(na, beta); % game-payoff matrix for each attack type by na and each sensor budge by beta
    for k = 1 : na % for each attack type
        As_sum = nchoosek(1 : N, alpha(k)); % all the scenarios for attack type k 
        psi = sdpvar(k, 1); % optimal variable for each attack type (corresponding to the number of attack nodes)
        BAnv = sdpvar(N, k); % attack entry matrix
        Qk = sdpvar(1, 1);
        F = [];
        F = [F, omegaM >= 0];
        F = [F, omegaM <= Vinf_hat * zMv];
        F = [F, psi >= 0];
        lmi = [-L'*P-P*L P*BAnv;
                BAnv'*P -diag(psi)] + diag([w.^2; zeros(k, 1)]) - diag([omegaM; zeros(k, 1)]);
        F = [F, lmi <= 0];
        Qk = delta' * omegaM + E * ones(1, k) * psi;
        R_cost = optimizer(F, Qk, options, {zMv, BAnv}, Qk);
        for As = 1 : size(As_sum, 1) % for each attack set of attack type k
            An = As_sum(As,:); % take an attack set
            BAn = [];
            for ka = 1:size(An,2)
                ea = zeros(N,1);
                ea(An(ka)) = 1;
                BAn = [BAn,ea]; % create attack entry matrix based on the taken attack set
            end 
            
            for kb = 1 : beta % for each monitor budget
                zMpos = nchoosek(1:N, kb); % all possible monitor sets
                for kr = 1 : size(zMpos,1) % for each monitor set
                    zMpick = zMpos(kr,:); % take a monitor set
                    zM = zeros(N,1); 
                    zM(zMpick) = 1; % create a monitor set based on the taken monitor set
                    % compute the worst-case disruption for a fixed pair of
                    % attack set and monitor set
                    tic
                    R_matrix{k, kb}(As, kr) = R_cost({zM, BAn}) + kappa' * zM; % fill in a cell of the game-payoff matrix 
                    exptimeBF(expk) = exptimeBF(expk) + toc; % computation time for this experiment
                end
            end
        end
    end

    % cost evaluation
    R_cost = cell(beta, 1);
    R_cost_min = zeros(beta, 1);
    for kb = 1 : beta
        R_cost{kb, 1} = zeros(size(max(R_matrix{1, kb})));
        for ka = 1 : na
            R_cost{kb, 1} = R_cost{kb, 1} + phi(ka) * max(R_matrix{ka, kb});      
        end
        expRd{expk} = [expRd{expk}, R_cost{kb, 1}]; % record all the worst-case disruptions by each admissible monitor set
        R_cost_min(kb) = min(R_cost{kb, 1}); % the minimum defense cost for each sensor budget  
    end
    plot(expk*ones(size(expRd{expk}, 2)-1,1),expRd{expk}(2:end),'o','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',4)
    hold on
end

%% Figure
set(gca,'Fontsize',16);
ylabel('Defense cost','FontSize',16);
xlabel('Experiment','FontSize',16);
xticks(0:5:expN);
xlim([0 expN+1]);
ylim([1 max(expRd)*1.3]);
grid;
legend('Optimal monitor set','All the monitor sets');
%print(gcf,'figs/10Nodeopt','-dpng','-r300');
