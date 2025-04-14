function [L,A] = ConnectedGraph(N,direct)
% ==========================================================
% ==========================================================
% This function generates a strongly connected graph
% Outputs are matrices L and A
% L is the Laplacian matrix reprents the graph with diagonal as in-degree
% matrix
% A is an adjacency matrix, a(i,j) if there exists an edge from j to i
% Inputs are N and direct
% N is the size of the network
% direct = 1 to generate a directed graph
% direct = 0 to generate an undirected graph
% ==========================================================
% ==========================================================
    L = zeros(N,N);
    A = zeros(N,N);
    edge_connect = [0, 0, 0, 1];
    while 1
        % Degree matrix or in-degree matrix
        D = zeros(N,N);
        for a = 1:1:N
            for m = 1:1:N
                D(a,a) = D(a,a) + A(a,m);
            end
        end
        % Laplacian matrix
        L = D - A;
        G = digraph(A');
        G_concomp = conncomp(G, 'Type', 'strong', 'OutputForm', 'cell');
        if size(G_concomp)*[0;1] == 1
            break
        end
        % Connect edges
        for i = 1:1:N-1
            for j = i+1:1:N
                if A(i,j) == 0
                    %A(i,j) = randi([0 1],1,1);
                    A(i,j) = edge_connect(randi([1 length(edge_connect)],1,1));
                    A(i,j) = A(i,j) * (1 + rand(1,1) * 2e-1);
                end
            end
        end
        if direct == 0
            for a = 1:1:N
                for m = 1:1:a
                    A(a,m) = A(m,a); 
                end
            end
        else
            for a = 1:1:N-1
                for m = i+1:1:N
                    if A(m,a) == 0
                        A(m,a) = edge_connect(randi([1 length(edge_connect)],1,1));
                        A(m,a) = A(m,a) * (1 + rand(1,1) * 2e-1);
                    end
                end
            end
        end

    end


    % while 1
    %     A = zeros(N,N);
    %     for i = 1:1:N-1
    %         for j = i+1:1:N
    %             A(i,j) = randi([0 1],1,1);
    %         end
    %     end
    %     %
    %     if direct == 0
    %         for a = 1:1:N
    %             for m = 1:1:a
    %                 A(a,m) = A(m,a); 
    %             end
    %         end
    %     else
    %         for a = 1:1:N-1
    %             for m = i+1:1:N
    %                 A(m,a) = randi([0,1],1,1);  
    %             end
    %         end
    %     end
    %     % Degree matrix or in-degree matrix
    %     D = zeros(N,N);
    %     for a = 1:1:N
    %         for m = 1:1:N
    %             D(a,a) = D(a,a) + A(a,m);
    %         end
    %     end
    %     % Laplacian matrix
    %     L = D - A;
    %     if direct == 0
    %         eL = eig(L);
    %         if eL(2) > 1e-3
    %             %disp('Graph is connected');
    %             break
    %         end
    %     else
    %         irA = eye(N);
    %         for i = 1:1:N-1
    %             irA = irA + (eye(N)-0.2*L)^i; 
    %         end
    %         if min(min(irA)) > 0
    %             %disp('Graph is connected');
    %             break
    %         end
    %     end
    % end
end